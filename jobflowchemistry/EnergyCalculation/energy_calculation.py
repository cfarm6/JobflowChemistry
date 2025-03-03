# from pydantic.dataclasses import dataclass
from dataclasses import dataclass, fields
from jobflow import Response, job, Maker

from ..utils import rdkit2ase, ase2rdkit
from rdkit.Chem import rdmolops, rdmolfiles, rdchem
from ..Calculators import ORCACalculator
from ..Calculators.ASECalculator import ASECalculator
from ..Calculators.TBLite import TBLiteCalculator
from ..Calculators.AimNet2 import AimNet2Calculator
from ..Calculators.xTBCalculator import xTBCalculator
from ..Calculators.PTBCalculator import PTBCalculator
from ..task_node import TaskNode
from typing import Union
from ..Structure import Structure
from ..outputs import Settings, Properties


@dataclass
class EnergyCalculation(Maker):
    name: str = "Single Point Calculation"

    def get_settings(self):
        raise NotImplementedError

    def calculate_energy(self, molecule: Structure):
        raise NotImplementedError

    @job(files="files", settings="settings", properties="properties")
    def make(self, structure: Structure):
        if type(structure) is list:
            jobs = [self.make(s) for s in structure]
            return Response(
                output={
                    "structure": [x.output["structure"] for x in jobs],
                    "settings": Settings({}),
                    "properties": [x.output["properties"] for x in jobs],
                },
                replace=jobs,
            )
        if structure.GetNumConformers() > 1:
            jobs = []
            for confId in range(structure.GetNumConformers()):
                s = Structure(rdchem.Mol(structure, confId=confId))
                jobs.append(self.make(s))
            return Response(
                output={
                    "structure": [x.output["structure"] for x in jobs],
                    "settings": Settings({}),
                    "properties": [x.output["properties"] for x in jobs],
                },
                replace=jobs,
            )

        properties = self.calculate_energy(structure)
        if "Global" in properties:
            for k, v in properties["Global"].items():
                if type(v) is list or type(v) is dict:
                    structure.SetProp(k, ",".join([str(x) for x in v]), computed=True)
                    continue
                structure.SetDoubleProp(k, float(v), computed=True)
        if "Atomic" in properties:
            for k, v in properties["Atomic"].items():
                for i, atom in enumerate(structure.GetAtoms()):
                    atom.SetDoubleProp(k, float(v[i]))
        if "Bond" in properties:
            for k, v in properties["Bond"].items():
                for i in v:
                    bond = structure.GetBondBetweenAtoms(i["atom1"], i["atom2"])
                    if bond is None:
                        continue
                    bond.SetDoubleProp(k, float(i["value"]))
        settings = self.get_settings()
        return Response(
            output={
                "structure": Structure(structure),
                "files": rdmolfiles.MolToV3KMolBlock(structure),
                "settings": Settings(settings),
                "properties": Properties(properties),
            },
            stored_data=properties,
        )


@dataclass
class ASEEnergyCalculator(ASECalculator, EnergyCalculation):
    name: str = "ASE Energy Calculator"

    def get_settings(self):
        raise NotImplementedError

    def calculate_energy(self, structure: Structure):
        atoms = self.set_calculator(structure)
        properties = self.get_properties(atoms)
        return properties


@dataclass
class TBLiteEnergyCalculation(TBLiteCalculator, ASEEnergyCalculator):
    name: str = "TBLite Energy Calculation"


@dataclass
class AimNet2EnergyCalculation(AimNet2Calculator, ASEEnergyCalculator):
    name: str = "AimNet2 Energy Calculation"


@dataclass
class ORCAEnergyCalculation(ORCACalculator, EnergyCalculation):
    name: str = "ORCA Energy Calculation"
    runtype: str = "ENERGY"

    def calculate_energy(self, structure: Structure):
        self.run_orca(structure)
        properties = self.get_properties(structure)
        return properties


@dataclass
class xTBEnergyCalculation(xTBCalculator, EnergyCalculation):
    name: str = "xTB Energy Calculation"

    def __post_init__(self):
        # List of prefixes to delete
        prefixes = ["md_", "opt_", "thermo_", "md_", "hess_", "modef_", "cube_"]

        # Remove attributes starting with any of the prefixes
        for field in fields(self):
            if any(field.name.startswith(prefix) for prefix in prefixes):
                delattr(self, field.name)

    def calculate_energy(self, structure: Structure):
        self.run_xtb(structure)
        properties = self.get_properties(structure)
        return properties

@dataclass
class PTBCalculation(PTBCalculator, EnergyCalculation):
    name: str = "PTB Energy Calculation"
    raman = False
    def calculate_energy(self, structure: Structure):
        properties = self.get_properties(structure)
        return properties
