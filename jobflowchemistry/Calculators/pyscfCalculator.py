from dataclasses import dataclass, field
from typing import Literal

# from dataclasses import dataclass, field
from .ASECalculator import ASECalculator
from rdkit.Chem import rdmolops, rdchem
from ..utils import rdkit2ase
from ase import Atoms
from pyscf import dft

xc_functionals = list(dft.XC.keys())

@dataclass
class PYSCFCalculator(ASECalculator):
    name: str = "PySCCF Calculator"
    xc_functional: str = None
    basiss_set: str = None
    grid: Literal[0, 1, 2, 3, 4, 5, 6, 7, 8, 9] = 3
    

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
            "atomic": {
                "Mulliken Partial Charges [e]": charges,
            },
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
