from pydantic.dataclasses import dataclass
from jobflow import Maker
from ase import Atoms
@dataclass
class ASECalculator():
    """
    Base class for ASE-based quantum chemistry calculators.

    Provides an interface for setting up, running, and extracting properties from ASE calculators.

    Attributes
    ----------
    name : str
        Name of the calculator.

    Methods
    -------
    set_calculator(atoms)
        Abstract method to set up the ASE calculator for the given atoms.
    get_properties(atoms)
        Abstract method to extract properties from the ASE Atoms object.
    get_settings()
        Abstract method to return calculator settings.
    """
    name: str = "ASE Calculator"
    def set_calculator(self, atoms: Atoms):
        raise NotImplementedError
    def get_properties(self, atoms: Atoms):
        raise NotImplementedError
    def get_settings(self):
        raise NotImplementedError