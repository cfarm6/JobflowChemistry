from pydantic.dataclasses import dataclass
from jobflow import Maker
from ase import Atoms
@dataclass
class ASECalculator():
    name: str = "ASE Calculator"
    def set_calculator(self, atoms: Atoms):
        raise NotImplementedError
    def get_properties(self, atoms: Atoms):
        raise NotImplementedError
    def get_settings(self):
        raise NotImplementedError