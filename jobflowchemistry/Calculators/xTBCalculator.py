from dataclasses import dataclass, field
from typing import Literal
from rdkit.Chem import rdmolfiles, rdmolops, rdchem
import json
import subprocess
import re


@dataclass
class xTBCalculator:
    name: str = "xTB Calculator"
    executable: str = "xtb"
    chrg: int = 0
    spin: int = 0
    cma: bool = True
    keyword_ceh: bool = False
    gfn_method: Literal[1, 2] = 2
    gfn_scc: bool = True
    gfn_periodic: bool = False
    gfn_dispscale: float = 1.0
    scc_maxiterations: int = 250
    scc_temp: float = 300.0
    scc_broydamp: float = 0.4
    opt_optlevel: Literal[
        "crude", "sloppy", "loose", "normal", "tight", "verytight", "extreme"
    ] = "normal"
    opt_engine: Literal["rf", "lbfgs", "inertial"] = "rf"
    opt_microcycle: int = 25
    opt_maxcycle: int = 0
    opt_hlow: float = 0.100e-1
    opt_s6: float = 20.0
    opt_kstretch: float = 0.4
    opt_kbend: float = 0.13
    opt_koutofp: float = 0.0
    opt_kvdw: float = 0.0
    opt_kes: float = 0.0
    opt_rcut: float = 8.3666002653407556
    opt_exact_rf: bool = False
    opt_average_conv: bool = False
    thermo_temp: float = 298.14999999999998
    thermo_sthr: float = 50.000000000000000
    thermo_imagthr: float = -20.000000000000000
    thermo_scale: float = 1.0000000000000000
    md_temp: float = 298.14999999999998
    md_time: float = 50.000000000000000
    md_dump: float = 50.000000000000000
    md_velo: int = 0
    md_nvt: Literal[1, 0] = 1
    md_skip: int = 500
    md_step: float = None
    md_hmass: int = None
    md_shake: int = None
    md_sccacc: float = 2.0000000000000000
    md_forcewrrestart: bool = False
    hess_sccacc: float = 0.29999999999999999
    hess_step: float = 0.50000000000000001e-2
    hess_scale: float = 1.0000000000000000
    modef_n: int = 31
    modef_step: float = 1.0000000000000000
    modef_updat: float = 0.20000000000000001
    modef_local: Literal[0, 1] = 0
    modef_vthr: float = 0.0000000000000000
    modef_prj: int = 0
    modef_mode: int = 0
    cube_step: float = 0.40000000000000002
    cube_pthr: float = 0.50000000000000003e-1
    cube_boff: float = 3.0
    solvent: Literal[
        "acetone",
        "acetonitrile",
        "aniline",
        "benzaldehyde",
        "benzene",
        "ch2cl2",
        "chcl3",
        "cs2",
        "dioxane",
        "dmf",
        "dmso",
        "ether",
        "ethylacetate",
        "furane",
        "hexandecane",
        "hexane",
        "methanol",
        "nitromethane",
        "octanol",
        "woctanol",
        "phenol",
        "toluene",
        "thf",
        "water",
    ] = None
    alpb: bool = False
    cosmo: bool = False
    cpcmx: Literal[
        "1,1,1-trichloroethane",
        "1,1,2-trichloroethane",
        "1,2,4-trimethylbenzene",
        "1,2-dibromoethane",
        "1,2-dichloroethane",
        "1,2-ethanediol",
        "1,4-dioxane",
        "1-bromo-2-methylpropane",
        "1-bromooctane",
        "1-bromopentane",
        "1-bromopropane",
        "1-butanol",
        "1-chlorohexane",
        "1-chloropentane",
        "1-chloropropane",
        "1-decanol",
        "1-fluorooctane",
        "1-heptanol",
        "1-hexanol",
        "1-hexene",
        "1-hexyne",
        "1-iodobutane",
        "1-iodohexadecane",
        "1-iodopentane",
        "1-iodopropane",
        "1-nitropropane",
        "1-nonanol",
        "1-octanol",
        "1-pentanol",
        "1-pentene",
        "1-propanol",
        "2,2,2-trifluoroethanol",
        "2,2,4-trimethylpentane",
        "2,4-dimethylpentane",
        "2,4-dimethylpyridine",
        "2,6-dimethylpyridine",
        "2-bromopropane",
        "2-butanol",
        "2-chlorobutane",
        "2-heptanone",
        "2-hexanone",
        "2-methoxyethanol",
        "2-methyl-1-propanol",
        "2-methyl-2-propanol",
        "2-methylpentane",
        "2-methylpyridine",
        "2-nitropropane",
        "2-octanone",
        "2-pentanone",
        "2-propanol",
        "2-propen-1-ol",
        "e-2-pentene",
        "3-methylpyridine",
        "3-pentanone",
        "4-heptanone",
        "4-methyl-2-pentanone",
        "4-methylpyridine",
        "5-nonanone",
        "acetic acid",
        "acetone",
        "acetonitrile",
        "acetophenone",
        "ammonia",
        "aniline",
        "anisole",
        "benzaldehyde",
        "benzene",
        "benzonitrile",
        "benzyl alcohol",
        "bromobenzene",
        "bromoethane",
        "bromoform",
        "butanal",
        "butanoic acid",
        "butanone",
        "butanonitrile",
        "butyl ethanoate",
        "butylamine",
        "n-butylbenzene",
        "sec-butylbenzene",
        "tert-butylbenzene",
        "carbon disulfide",
        "carbon tetrachloride",
        "chlorobenzene",
        "chloroform",
        "a-chlorotoluene",
        "o-chlorotoluene",
        "conductor",
        "m-cresol",
        "o-cresol",
        "cyclohexane",
        "cyclohexanone",
        "cyclopentane",
        "cyclopentanol",
        "cyclopentanone",
        "decalin",
        "cis-decalin",
        "n-decane",
        "dibromomethane",
        "dibutylether",
        "o-dichlorobenzene",
        "e-1,2-dichloroethene",
        "z-1,2-dichloroethene",
        "dichloromethane",
        "diethyl ether",
        "diethyl sulfide",
        "diethylamine",
        "diiodomethane",
        "diisopropyl ether",
        "cis-1,2-dimethylcyclohexane",
        "dimethyl disulfide",
        "n,n-dimethylacetamide",
        "n,n-dimethylformamide",
        "dimethylsulfoxide",
        "diphenylether",
        "dipropylamine",
        "n-dodecane",
        "ethanethiol",
        "ethanol",
        "ethyl acetate",
        "ethyl methanoate",
        "ethyl phenyl ether",
        "ethylbenzene",
        "fluorobenzene",
        "formamide",
        "formic acid",
        "furan",
        "n-heptane",
        "n-hexadecane",
        "n-hexane",
        "hexanoic acid",
        "iodobenzene",
        "iodoethane",
        "iodomethane",
        "isopropylbenzene",
        "p-isopropyltoluene",
        "mesitylene",
        "methanol",
        "methyl benzoate",
        "methyl butanoate",
        "methyl ethanoate",
        "methyl methanoate",
        "methyl propanoate",
        "n-methylaniline",
        "methylcyclohexane",
        "n-methylformamide",
        "nitrobenzene",
        "nitroethane",
        "nitromethane",
        "o-nitrotoluene",
        "n-nonane",
        "n-octane",
        "n-pentadecane",
        "octanol(wet)",
        "pentanal",
        "n-pentane",
        "pentanoic acid",
        "pentyl ethanoate",
        "pentylamine",
        "perfluorobenzene",
        "phenol",
        "propanal",
        "propanoic acid",
        "propanonitrile",
        "propyl ethanoate",
        "propylamine",
        "pyridine",
        "tetrachloroethene",
        "tetrahydrofuran",
        "tetrahydrothiophene-s,s-dioxide",
        "tetralin",
        "thiophene",
        "thiophenol",
        "toluene",
        "trans-decalin",
        "tributylphosphate",
        "trichloroethene",
        "triethylamine",
        "n-undecane",
        "water",
        "xylene",
        "m-xylene",
        "o-xylene",
        "p-xylene",
        "",
    ] = None
    write_esp: bool = False
    write_mos: bool = False
    write_gbw: bool = False
    write_tm_mos: bool = False
    write_tm_basis: bool = False
    write_lmo: bool = False
    write_density: bool = False
    write_spin_population: bool = False
    write_spin_density: bool = False
    write_fod: bool = False
    write_fod_population: bool = False
    write_wiberg: bool = True
    write_wbo_fragments: bool = False
    write_dipole: bool = True
    write_charges: bool = True
    write_mulliken: bool = True
    write_orbital_energies: bool = True
    write_inertia: bool = True
    write_distances: bool = True
    write_angles: bool = False
    write_torsions: bool = False
    write_final_struct: bool = True
    write_geosum: bool = True
    write_stm: bool = False
    write_modef: bool = False
    write_gbsa: bool = False
    write_json: bool = True

    # Theory
    def __post_init__(self):
        if self.solvent and self.cpcmx:
            raise ValueError(
                "Both CPCMX and ALPB/COSMO Solvation models cannot be used simultaneously"
            )

    def get_settings(self):
        settings = {}
        for key in self.__dict__.keys():
            settings[key] = getattr(self, key)
        return settings

    def get_properties(self, molecule: rdchem.Mol):
        properties = {"Global": {}, "Atomic": {}, "Bond": {}}
        with open("xtbout.json", "r") as f:
            _props = json.load(f)
        properties["Global"]["HOMO-LUMO Gap [eV]"] = _props["HOMO-LUMO gap / eV"]
        properties["Global"]["Total Energy [eV]"] = _props["total energy"]
        properties["Global"]["Electronic Energy [eV]"] = _props["electronic energy"]
        properties["Global"]["Dipole [au]"] = {
            "x": _props["dipole / a.u."][0],
            "y": _props["dipole / a.u."][1],
            "z": _props["dipole / a.u."][2],
        }
        properties["Global"]["Orbitals"] = [
            {"Occupation [-]": x[1], "Orbital Energy [eV]": x[0]}
            for x in zip(
                _props["orbital energies / eV"], _props["fractional occupation"]
            )
        ]
        properties["Atomic"]["Mulliken Partial Charges [e]"] = _props["partial charges"]
        properties["Bond"]["Wiberg Bond Order [-]"] = []
        if self.keyword_ceh:
            properties["Atomic"]["Charge Extended Huckel Charge [e]"] = []
            with open("ceh.charges") as f:
                for line in f.readlines():
                    if re.match(r"^\s*\d", line):
                        properties["Atomic"][
                            "Charge Extended Huckel Charge [e]"
                        ].append(float(line.split()[0]))
        with open("wbo", "r") as f:
            wbo = f.readlines()
            for line in wbo:
                atom1, atom2, _wbo = line.split()
                properties["Bond"]["Wiberg Bond Order [-]"].append(
                    {
                        "atom1": int(atom1) - 1,
                        "atom2": int(atom2) - 1,
                        "value": float(_wbo),
                    }
                )
        return properties

    def set_keywords(self):
        self.keywords = [self.executable, "input.mol", "--input", "xtb.input", "--molden"]
        return

    def run_xtb(self, molecule: rdchem.Mol):
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

        self.set_keywords()
        # create the xtb input file
        with open("xtb.input", "w") as f:
            f.write(f"$chrg {self.chrg}\n")
            f.write(f"$spin {self.spin}\n")
            if self.cma:
                f.write(f"$cma\n")
            grouped_fields = {}

            for attr, value in self.__dict__.items():
                if "_" in attr:  # Only consider fields with "_"
                    key, subkey = attr.split("_", 1)
                    if key == "keyword":
                        continue
                    if key == "solvation":
                        if self.solvation_solvent is None:
                            continue
                    if key not in grouped_fields:
                        grouped_fields[key] = []
                    grouped_fields[key].append(f"{subkey.replace("_", " ")}={value}")
            document_lines = []
            for key, subkeys in grouped_fields.items():
                document_lines.append(f"${key}")
                document_lines.extend(subkeys)
            f.write("\n".join(document_lines))
            f.write("\n$end")

        if self.solvent is not None:
            if self.cosmo:
                self.keywords.append(f"--cosmo {self.solvent}")
            elif self.alpb:
                self.keywords.append(f"--alpb {self.solvent}")
            else:
                raise ValueError("ALPB or COSMO must be specified")
        if self.cpcmx:
            self.keywords.append(f"--cpcmx {self.cpcmx}")
        if self.keyword_ceh:
            self.keywords.append("--ceh")
        command = " ".join(self.keywords)
        subprocess.call(f"{command} > log.out", shell=True)
        return None
