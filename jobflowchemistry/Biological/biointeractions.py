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
class BoltzFolding(Maker):
    name: str = "Boltzfolding"

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
