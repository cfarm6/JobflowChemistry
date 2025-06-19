from pydantic.dataclasses import dataclass

from .ASECalculator import ASECalculator
from rdkit.Chem import rdmolops
from ..utils import rdkit2ase
from ase import Atoms
from rdkit.Chem import rdchem


@dataclass
class PETMADCalculator(ASECalculator):
    name: str = "PET-MAD Calculator"
    version: str = "latest"
    device: str = "cuda"

    def get_settings(self):
        settings = {"version": self.model, "device": self.device}

    def get_properties(self, atoms: Atoms):
        energy = atoms.get_total_energy()
        forces = atoms.get_forces()
        properties = {
            "Global": {"Total Energy [eV]": energy},
            "Atomic": {"Forces [eV/A]": forces},
        }
        return properties

    def set_calculator(self, molecule: rdchem.Mol):
        from pet_mad.calculator import PETMADCalculator as PMCalc

        atoms = rdkit2ase(molecule)
        atoms.calc = PMCalc(version=self.version, device=self.device)
        Warning(
            "WARNING: This model does not account for molecular charge and multiplicity"
        )
        return atoms
