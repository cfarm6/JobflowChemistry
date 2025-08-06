# from pydantic import Field
from dataclasses import dataclass, field

# from pydantic.dataclasses import dataclass

from jobflow import Response, job, Flow, Maker
import pickle
from rdkit.Chem import rdmolops, rdmolfiles, rdchem
from rdkit import Chem
from typing import Any, Tuple, Union, Dict
import os
from typing import Literal
import subprocess

from ..Calculators.PTBCalculator import PTBCalculator
from ..Structure import Structure
from ..outputs import Properties, Settings


@dataclass
class SpectraCalculator(Maker):
    """
    Base class for spectra calculation jobs in workflows.

    Provides an interface for calculating molecular spectra (e.g., IR, Raman).

    Attributes
    ----------
    name : str
        Name of the calculator.

    Methods
    -------
    get_settings()
        Abstract method to return calculator settings.
    get_properties(molecule)
        Abstract method to compute spectra properties for a molecule.
    make(structure)
        Job for executing the spectra calculation and returning results as a Response.
    """
    name: str = "Spectra Calculation"

    def get_settings(self):
        raise NotImplementedError

    def get_properties(self, molecule: Structure):
        raise NotImplementedError

    @job(files="files", settings="settings", properties="properties")
    def make(self, structure: Structure):
        # molecule = pickle.loads(molecule)
        if type(structure) is list:
            jobs = [self.make(s) for s in structure]
            return Response(
                output={
                    "structure": [x.output["structure"] for x in jobs],
                    "settings": Settings({}),
                    "files": [x.output["files"] for x in jobs],
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
                    "files": [x.output["files"] for x in jobs],
                    "settings": Settings({}),
                    "properties": [x.output["properties"] for x in jobs],
                },
                replace=jobs,
            )
        properties = self.get_properties(structure)

        if "global" in properties:
            for k, v in properties["global"].items():
                if type(v) is list: continue
                structure.SetDoubleProp(k, float(v), computed=True)
        if "atomic" in properties:
            for k, v in properties["atomic"].items():
                for i, atom in enumerate(structure.GetAtoms()):
                    if type(v[i]) is list:
                        continue
                    atom.SetDoubleProp(k, float(v[i]))
        if "bond" in properties:
            for k, v in properties["bond"].items():
                for i in v:
                    bond = structure.GetBondBetweenAtoms(i["atom1"], i["atom2"])
                    if bond is None:
                        continue
                    bond.SetDoubleProp(k, float(i["value"]))
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


@dataclass
class PTBRamanSpectra(PTBCalculator, SpectraCalculator):
    """
    Job for calculating Raman spectra using the PTB (xtb) calculator.

    Attributes
    ----------
    name : str
        Name of the job.
    raman : bool
        Whether to compute Raman spectra (always True for this class).
    """
    name: str = "PTB Raman Spectra"
    raman: bool = True
    pass
