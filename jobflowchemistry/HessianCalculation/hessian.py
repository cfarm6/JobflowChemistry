from ..utils import ase2rdkit
from ..Calculators import ORCACalculator
from ..Calculators.TBLite import TBLiteCalculator
from ..Calculators.ASECalculator import ASECalculator
from ..Calculators.AimNet2 import AimNet2Calculator
from ..Calculators.xTBCalculator import xTBCalculator
from ..Structure import Structure
from ..outputs import Settings, Properties
from ..utils import parse_gaussian_output

from dataclasses import field, fields
from dataclasses import dataclass

from jobflow import job, Response, Maker
from rdkit.Chem import rdmolfiles, rdchem

@dataclass
class HessianCalculation(Maker):
    name: str = "Hessian Calculation"

    def calculate_hessian(self, structure: Structure):
        raise NotImplementedError

    def get_settings(self):
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

        properties = self.calculate_hessian(structure)
        settings = self.get_settings()

        if "global" in properties:
            for k, v in properties["global"].items():
                if type(v) is list:
                    continue
                structure.SetDoubleProp(k, float(v), computed=True)
        if "atomic" in properties:
            for k, v in properties["atomic"].items():
                for i, atom in enumerate(structure.GetAtoms()):
                    atom.SetDoubleProp(k, float(v[i]))
        if "bond" in properties:
            for k, v in properties["bond"].items():
                for i in v:
                    bond = structure.GetBondBetweenAtoms(i["atom1"], i["atom2"])
                    if bond is None:
                        continue
                    bond.SetDoubleProp(k, float(i["value"]))
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
class xTBHessian(xTBCalculator, HessianCalculation):
    name: str = "xTB Hessian"
    preoptimize: bool = False
    def __post_init__(self):
        # List of prefixes to delete
        prefixes = ["md_", "thermo_", "md_", "hess_", "modef_", "cube_"]
        if not self.preoptimize:
            prefixes.append("--opt")

        # Remove attributes starting with any of the prefixes
        for field in fields(self):
            if any(field.name.startswith(prefix) for prefix in prefixes):
                delattr(self, field.name)

    def set_keywords(self):
        super().set_keywords()
        if self.preoptimize:
            self.keywords.append("--ohess")
        else:
            self.keywords.append("--hess")
        return

    def get_properties(self, molecule):
        properties = super().get_properties(molecule)
        (
            frequencies,
            reduced_masses,
            force_constants,
            ir_intensities,
            raman_activity,
            depolarization,
            displacement_vectors,
        ) = parse_gaussian_output("g98.out")
        properties["spectra"] = {"Normal Modes": {
            "Frequency [cm^-1]": frequencies,
            "Reduced Masses [amu]": reduced_masses,
            "Force Constants [mDyne/A]": force_constants,
            "IR Intensities [km/mol]": ir_intensities,
            "Raman Activity [A^4/amu]": raman_activity,
            "Displacements": displacement_vectors
        }}
        return properties

    def calculate_hessian(self, structure: Structure):
        self.set_keywords()
        self.run_xtb(structure)
        properties = self.get_properties(structure)
        return properties
