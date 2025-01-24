from ..utils import ase2rdkit
from ..Calculators import ORCACalculator
from ..Calculators.TBLite import TBLiteCalculator
from ..Calculators.ASECalculator import ASECalculator
from ..Calculators.AimNet2 import AimNet2Calculator
from ..Calculators.xTBCalculator import xTBCalculator

from dataclasses import field, fields
from dataclasses import dataclass
# from pydantic.dataclasses import dataclass
from typing import Literal

from jobflow import job, Response, Maker
from rdkit.Chem import rdmolfiles, rdchem
import ase.optimize

from ..Structure import Structure
from ..outputs import Settings, Properties

@dataclass
class GeometryOptimization(Maker):
    name: str = "Geometry Optimization"

    def optimize_structure(self, structure: Structure):
        raise NotImplementedError

    @job(files="files", settings="settings", properties="properties")
    def make(self, structure: Structure):
        if type(structure) is list:
            jobs = [self.make(s) for s in structure]
            return Response(
                output={
                    "structure": [x.output["structure"] for x in jobs],
                    "settings": Settings({}),
                    "properties": [x.output["properties"] for x in jobs],
                },
                replace=jobs,
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
                    "properties": [x.output["properties"] for x in jobs],
                },
                replace=jobs,
            )
        structure, properties, settings = self.optimize_structure(structure)
        if "global" in properties:
            for k,v in properties["global"].items():
                if type(v) is list: continue
                structure.SetDoubleProp(k, float(v), computed=True)
        if "atomic" in properties:
            for k,v in properties["atomic"].items():
                for i, atom in enumerate(structure.GetAtoms()):
                    atom.SetDoubleProp(k, float(v[i]))
        if "bond" in properties:
            for k,v in properties["bond"].items():
                for i in v:
                    bond = structure.GetBondBetweenAtoms(i['atom1'], i['atom2'])
                    if bond is None: continue
                    bond.SetDoubleProp(k, float(i['value']))

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


@dataclass
class xTBOptimization(xTBCalculator, GeometryOptimization):
    name: str = "xTB Optimizer"
    def __post_init__(self):
        # List of prefixes to delete
        prefixes = ["md_", "thermo_", "md_", "hess_", "modef_", "cube_"]

        # Remove attributes starting with any of the prefixes
        for field in fields(self):
            if any(field.name.startswith(prefix) for prefix in prefixes):
                delattr(self, field.name)
    def set_keywords(self):
        super().set_keywords()
        self.keywords.append("--opt")
        return
    def optimize_structure(self, structure: Structure):
        self.set_keywords()
        self.run_xtb(structure)
        properties = self.get_properties(structure)
        mol = rdmolfiles.MolFromMolFile("xtbopt.mol", sanitize=False, removeHs=False)
        settings = self.get_settings()
        return mol, properties, settings
