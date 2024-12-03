from ..utils import ase2rdkit
from ..Calculators import ORCACalculator
from ..Calculators.TBLite import TBLiteCalculator
from ..Calculators.ASECalculator import ASECalculator
from ..Calculators.AimNet2 import AimNet2Calculator

from dataclasses import field
from dataclasses import dataclass
# from pydantic.dataclasses import dataclass
from typing import Literal

from jobflow import job, Response, Maker
from rdkit.Chem import rdmolfiles, rdchem
import ase.optimize
import pickle
import subprocess
from icecream import ic

from ..Structure import Structure
from ..outputs import Settings, Properties

@dataclass
class GeometryOptimization(Maker):
    name: str = "Geometry Optimization"
    
    def optimize_structure(self, structure: Structure):
        raise NotImplementedError

    @job(files="files", settings="settings", properties="properties")
    def make(self, structure: Structure):
        structure, properties, settings = self.optimize_structure(structure)
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
class ASEOptimization(GeometryOptimization, ASECalculator):
    name: str = "ASE Optimization"
    optimizer: Literal["LBFGS", "BFGS", "FIRE", "FIRE2"] = field(default="LBFGS")
    fmax: float = 0.05
    steps: int = 250000

    def optimize_structure(self, structure: rdchem.Mol):
        atoms = self.set_calculator(structure)
        opt = getattr(ase.optimize, self.optimizer)(atoms, logfile="log.out")
        opt.run(self.fmax, self.steps)
        structure = ase2rdkit(atoms, structure)
        properties = self.get_properties(atoms)
        settings = self.get_settings()
        settings["fmax"] = self.fmax
        settings["steps"] = self.steps
        return structure, properties, settings


@dataclass
class TBLiteOptimization(ASEOptimization, TBLiteCalculator):
    name: str = "TBLite Optimization"


@dataclass
class AimNet2Optimization(ASEOptimization, AimNet2Calculator):
    name: str = "AimNet2 Optimization"

@dataclass
class ORCAOptimization(GeometryOptimization, ORCACalculator):
    name: str = "ORCA Optimization"
    runtype: Literal["OPT"] = "OPT"

    def optimize_structure(self, molecule: rdchem.Mol):
        super().run_orca(molecule)
        output_structure_filename = f"{self.base_name}.xyz"
        atoms = ase.io.read(output_structure_filename)
        molecule = ase2rdkit(atoms, molecule)
        settings = super().get_settings()
        properties = super().get_properties(molecule)
        return molecule, properties, settings
