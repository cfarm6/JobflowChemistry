from dataclasses import dataclass, field
from typing import Literal
from rdkit.Chem import rdmolfiles, rdmolops, rdchem
import json
import subprocess
import re
import os
from ..utils import fix_json_file
@dataclass
class PTBCalculator():
    name: str = "PTB Calculator"
    executable: str = "xtb"
    chrg: int = 0
    spin: int = 0
    cma: bool = True
    ceh: bool = False
    raman: bool = False

    # Theory
    def get_settings(self):
        d = {}
        for key, value in vars(self).items():
            if key == "name" or value is None:
                continue
            d[key] = value
        return d

    def get_properties(self, molecule: rdchem.Mol):
        self._run_ptb(molecule)
        properties = {"Global": {}, "Atomic": {}, "Bond": {}}
        fix_json_file("xtbout.json", "xtbout.json")
        with open("xtbout.json", "r") as f:
            _props = json.load(f)
        properties["Global"]["Total Energy [eV]"] = _props["total energy"]
        properties["Global"]["HOMO-LUMO Gap [eV]"] = _props["HOMO-LUMO gap / eV"]
        properties["Global"]["Electronic Energy [eV]"] = _props["electronic energy"]
        properties["Global"]["Dipole [au]"] = _props["dipole / a.u."]
        properties["Global"]["Orbital Energies [eV]"] = _props["orbital energies / eV"]
        properties["Global"]["Occupations [-]"] = _props["fractional occupation"]
        properties["Atomic"]["Partial Charges [e]"] = _props["partial charges"]
        if self.ceh:
            properties["Atomic"]["Charge Extended Huckel Charge [e]"] = []
            with open("ceh.charges") as f:
                for line in f.readlines():
                    if re.match(r"^\s*\d", line):
                        properties["Atomic"][
                            "Charge Extended Huckel Charge [e]"
                        ].append(float(line.split()[0]))
        if "bond orders" in _props:
            properties["Bond"]["Bond Order [-]"] = []
            for bond in _props["bond orders"]:
                atom1, atom2, _wbo = bond
                properties["Bond"]["Bond Order [-]"].append(
                    {
                        "atom1": int(atom1) - 1,
                        "atom2": int(atom2) - 1,
                        "value": float(_wbo),
                    }
                )
        if os.path.exists("vibspectrum"):
            properties["spectra"] = {}
            properties["spectra"]["Frequency [cm^-1]"] = _props[ "vibrational frequencies / rcm"]
            properties["spectra"]["Reduced Masses [amu]"] = _props["reduced masses"]
            properties["spectra"]["IR Intensity [km/mol]"] = _props[
                "IR intensities / km/mol"
            ]
            properties["spectra"]["Raman Activity [A^4/amu]"] = _props[
                "Raman activities / A^4/amu"
            ]
        return properties

    def _run_ptb(self, molecule: rdchem.Mol):
        keywords = [self.executable, "input.mol", "--ptb"]
        rdmolfiles.MolToV3KMolFile(molecule, "input.mol")
        if self.chrg is None:
            self.chrg = rdmolops.GetFormalCharge(molecule)
        if self.spin is None:
            total_electrons = 0
            if self.chrg is not None:
                total_electrons += abs(self.chrg)
            for atom in molecule.GetAtoms():
                total_electrons += atom.GetNumRadicalElectrons()
            total_electron_spin = total_electrons / 2
            spin = int(2 * total_electron_spin + 1)
            if spin == 0:
                spin = 2
            self.spin = spin
        if self.ceh:
            keywords.append("--ceh")
        if self.raman:
            keywords.append("--raman")
        if self.chrg:
            keywords.append(f"--chrg {self.chrg}")
        if self.spin:
            keywords.append(f"--spin {self.spin}")
        keywords.append("--json")
        command = " ".join(keywords)
        subprocess.call(f"{command} > log.out", shell=True)
        return None
