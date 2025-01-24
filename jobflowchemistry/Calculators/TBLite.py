from dataclasses import dataclass, field
from typing import Literal

# from dataclasses import dataclass, field
from .ASECalculator import ASECalculator
from rdkit.Chem import rdmolops, rdchem
from ..utils import rdkit2ase
from ase import Atoms

@dataclass
class TBLiteCalculator(ASECalculator):
    name: str = "TBLite Calculator"
    charge: int = None
    multiplicity: int = None
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = field(default="GFN2-xTB")
    accuracy: float = 1.0
    electronic_temperature: float = 300.0
    max_iterations: int = 250
    initial_guess: Literal["sad", "eeq"] = field(default="sad")
    mixer_damping: float = 0.4
    electric_field: list = None
    spin_polarization: float = None
    verbosity: int = 0

    def get_settings(self):
        settings = {
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "method": self.method,
            "accuracy": 1.0,
            "electronic_temperature": self.electronic_temperature,
            "max_iterations": self.max_iterations,
            "initial_guess": self.initial_guess,
            "mixer_damping": self.mixer_damping,
            "electric_field": self.electric_field,
            "spin_polarization": self.spin_polarization,
        }
        return settings

    def get_properties(self, atoms: Atoms):
        total_energy = atoms.get_total_energy()
        charges = atoms.get_charges()
        dipole = atoms.get_dipole_moment()
        properties = {
            "global": {
                "Total Energy [eV]": total_energy,
                "Dipole Moment [D]": list(dipole),
            },
            "atomic":{
                "Mulliken Partial Charges [e]": charges,
            }
        }
        return properties

    def set_calculator(self, molecule: rdchem.Mol):
        from tblite.ase import TBLite

        if self.charge is None:
            self.charge = rdmolops.GetFormalCharge(molecule)
        atoms = rdkit2ase(molecule)
        atoms.calc = TBLite(
            method=self.method,
            charge=self.charge,
            multiplicity=self.multiplicity,
            accuracy=self.accuracy,
            electronic_temperature=self.electronic_temperature,
            max_iterations=self.max_iterations,
            initial_guess=self.initial_guess,
            mixer_damping=self.mixer_damping,
            electric_field=self.electric_field,
            spin_polarization=self.spin_polarization,
            verbosity=self.verbosity,
        )
        return atoms
