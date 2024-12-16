# from pydantic.dataclasses import dataclass
from dataclasses import dataclass
from jobflow import Response, job, Maker
import pickle
import ase
from icecream import ic
from ..utils import rdkit2ase, ase2rdkit
from rdkit.Chem import rdmolops, rdmolfiles, rdchem
import subprocess
from ..Calculators import ORCACalculator
from ..Calculators.ASECalculator import ASECalculator
from ..Calculators.TBLite import TBLiteCalculator
from ..Calculators.AimNet2 import AimNet2Calculator
from ..task_node import TaskNode
from typing import Union
from ..Structure import Structure
from ..outputs import Settings, Properties


@dataclass
class EnergyCalculation(Maker):
    name: str = "Energy Calculator"

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
                    "files": [x.output["files"] for x in jobs],
                    "properties": [x.output["properties"] for x in jobs],
                },
                addition=jobs,
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
                    "files": [x.output["files"] for x in jobs],
                    "properties": [x.output["properties"] for x in jobs],
                },
                addition=jobs,
            )

        properties = self.calculate_energy(structure)
        if "global" in properties:
            for k,v in properties["global"].items():
                if type(v) is list: 
                    structure.SetProp(k, ",".join([str(x) for x in v]), computed=True)
                    continue
                structure.SetDoubleProp(k, float(v), computed=True)
        if "atomic" in properties:
            for k,v in properties["atomic"].items():
                for i, atom in enumerate(structure.GetAtoms()):
                    atom.SetDoubleProp(k, float(v[i]))
        if "bond" in properties:
            for k,v in properties["bond"].items():
                bond = structure.GetBondBetweenAtoms(v[0], v[1])
                bond.SetDoubleProp(k, float(v[2]))
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
