# from dataclasses import dataclass, field
from pydantic.dataclasses import dataclass

from .ASECalculator import ASECalculator
from rdkit.Chem import rdmolops
from ..utils import rdkit2ase
from ase import Atoms
from rdkit.Chem import rdchem

@dataclass
class AimNet2Calculator(ASECalculator):
    name: str = "AimNet2 Calculator"
    charge: int = None
    multiplicity: int = None
    model: str = "aimnet2"

    def get_settings(self):
        settings = {
            "model": "aimnet2",
            "charge": self.charge,
            "multiplicity": self.multiplicity,
        }
        return settings

    def get_properties(self, atoms: Atoms):
        energy = atoms.get_total_energy()[0]
        charge = atoms.get_charges()
        properties = {
            "global": {"Total Energy [eV]": energy},
            "atomic": {"AimNet2 Partial Charges [e]": charge},
        }

        return properties

    def set_calculator(self, molecule: rdchem.Mol):
        from aimnet2calc import AIMNet2ASE

        if self.charge is None:
            self.charge = rdmolops.GetFormalCharge(molecule)
        atoms = rdkit2ase(molecule)
        atoms.calc = AIMNet2ASE(self.model)
        if self.charge is not None:
            atoms.calc.set_charge(self.charge)
        if self.multiplicity is not None:
            atoms.calc.set_mult(self.multiplicity)
        return atoms
