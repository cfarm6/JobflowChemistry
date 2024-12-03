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
        properties = self.calculate_energy(structure)
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
class ASEEnergyCalculator(
    ASECalculator, EnergyCalculation
):
    name: str = "ASE Energy Calculator"
    def get_settings(self):
        raise NotImplementedError

    def calculate_energy(self, structure: Structure):
        atoms = self.set_calculator(structure)
        properties = self.get_properties(atoms)
        return properties


@dataclass
class TBLiteEnergyCalculation(
    TBLiteCalculator, ASEEnergyCalculator
):
    name: str = "TBLite Energy Calculation"


@dataclass
class AimNet2EnergyCalculation(
    AimNet2Calculator, ASEEnergyCalculator
):
    name: str = "AimNet2 Energy Calculation"


@dataclass
class ORCAEnergyCalculation(
    ORCACalculator, EnergyCalculation
):
    name: str = "ORCA Energy Calculation"
    runtype: str = "ENERGY"

    def calculate_energy(self, structure: Structure):
        self.run_orca(structure)
        properties = self.get_properties(structure)
        return properties
