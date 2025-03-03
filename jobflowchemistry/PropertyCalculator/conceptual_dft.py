from ..outputs import Properties, Settings
from ..Structure import Structure
from .property_calculator import PropertyCalculator

from dataclasses import dataclass, field
from jobflow import Response, job, Flow, Maker
from rdkit.Chem import rdmolfiles, rdmolops
import subprocess
import re


@dataclass
class ConceptualDFTProperties(PropertyCalculator):
    name: str = "Conceptual DFT Properties"
    pass


@dataclass
class CalculateGlobalConceptualDFTProperties(ConceptualDFTProperties):
    name: str = "Calculate Global CDFT Properties"

    @job(files="files", settings="settings", properties="properties")
    def make(
        self,
        properties_neutral: Properties,
        properties_cation: Properties,
        properties_anion: Properties,
    ):
        neutral_state_energy = properties_neutral["Global"]["Total Energy [eV]"]
        cation_state_energy = properties_cation["Global"]["Total Energy [eV]"]
        anion_state_energy = properties_anion["Global"]["Total Energy [eV]"]
        vertical_ionization_potential = anion_state_energy - neutral_state_energy
        vertical_eletron_affinity = neutral_state_energy - cation_state_energy
        mulliken_electronegativity = (
            vertical_ionization_potential + vertical_eletron_affinity
        ) / 2
        chemical_potential = -mulliken_electronegativity
        hardness = vertical_ionization_potential - vertical_eletron_affinity
        softness = 1 / hardness
        elctrophilicity_index = chemical_potential**2 / 2 / hardness

        properties = {
            "Global": {
                "Vertical Ionization Potential [eV]": vertical_ionization_potential,
                "Vertical Electron Affinity [eV]": vertical_eletron_affinity,
                "Mulliken Electronegativity [eV]": mulliken_electronegativity,
                "Chemical Potential [eV]": chemical_potential,
                "Hardness [eV]": hardness,
                "Softness [1/eV]": softness,
                "Electrophilicity Index [eV]": elctrophilicity_index,
            }
        }
        return Response(
            output={
                "properties": Properties(properties),
            },
            stored_data=properties,
        )


@dataclass
class xTBVerticalElectronAffinityAndIonizationPotential(ConceptualDFTProperties):
    name: str = "xTB Vertical EA and IP"
    charge: int = None

    def get_settings(self):
        return {}

    def get_properties(self, structure: Structure):
        # Write Mol file
        rdmolfiles.MolToMolFile(structure, "input.mol")
        # Get Formal Charge
        if self.charge is None:
            self.charge = rdmolops.GetFormalCharge(structure)
        # Create the command
        command = [
            "xtb",
            "input.mol",
            "--vipea",
            "--chrg",
            str(self.charge),
            ">",
            "log.out",
        ]
        # Execute!
        subprocess.call(" ".join(command), shell=True)
        # Process output file
        properties = {}
        with open("log.out", "r") as f:
            for line in f:
                if line.startswith("delta SCC EA (eV)") or line.startswith(
                    "delta SCC IP (eV)"
                ):
                    results = line.split(":")
                    results[1] = float(results[1])
                    results[0] = (
                        results[0]
                        .replace("delta SCC", "vertical")
                        .replace("(", "[")
                        .replace(")", "]")
                        .title()
                        .replace("Ea", "Electron Affinity")
                        .replace("Ip", "Ionization Potential")
                        .replace("Ev", "eV")
                    )
                    properties[results[0]] = results[1]
        vertical_ionization_potential = properties["Vertical Electron Affinity [eV]"]
        vertical_eletron_affinity = properties["Vertical Ionization Potential [eV]"]
        mulliken_electronegativity = (
            vertical_ionization_potential + vertical_eletron_affinity
        ) / 2
        chemical_potential = -mulliken_electronegativity
        hardness = vertical_ionization_potential - vertical_eletron_affinity
        softness = 1 / hardness
        elctrophilicity_index = chemical_potential**2 / 2 / hardness

        return structure, {
            "Global": {
                "Vertical Ionization Potential [eV]": vertical_ionization_potential,
                "Vertical Electron Affinity [eV]": vertical_eletron_affinity,
                "Mulliken Electronegativity [eV]": mulliken_electronegativity,
                "Chemical Potential [eV]": chemical_potential,
                "Hardness [eV]": hardness,
                "Softness [1/eV]": softness,
                "Electrophilicity Index [eV]": elctrophilicity_index,
            }
        }


@dataclass
class xTBFukuiIndices(ConceptualDFTProperties):
    name: str = "xTB Fukui Indices"
    charge: int = None

    def get_settings(self):
        return {}

    def get_properties(self, structure: Structure):
        # Write Mol file
        rdmolfiles.MolToMolFile(structure, "input.mol")
        # Get Formal Charge
        if self.charge is None:
            self.charge = rdmolops.GetFormalCharge(structure)
        # Create the command
        command = [
            "xtb",
            "input.mol",
            "--vfukui",
            "--chrg",
            str(self.charge),
            ">",
            "log.out",
        ]
        # Execute!
        subprocess.call(" ".join(command), shell=True)
        # Process output file
        properties = {"Atomic": {}}
        in_fukui = False
        properties["Atomic"]["Fukui + [-]"] = []
        properties["Atomic"]["Fukui - [-]"] = []
        properties["Atomic"]["Fukui 0 [-]"] = []
        pattern = re.compile(r"^\s*(\d+\w)\s+([\d\.-]+)\s+([\d\.-]+)\s+([\d\.-]+)")
        with open("log.out", "r") as f:
            for line in f:
                if line.lstrip().startswith("#        f(+)     f(-)     f(0)"):
                    in_fukui = True
                    continue
                if in_fukui:
                    match = pattern.match(line)
                    if match:
                        properties["Atomic"]["Fukui + [-]"].append(
                            float(match.group(2))
                        )
                        properties["Atomic"]["Fukui - [-]"].append(
                            float(match.group(3))
                        )
                        properties["Atomic"]["Fukui 0 [-]"].append(
                            float(match.group(4))
                        )
                    else:
                        in_fukui = False
        return structure, properties
