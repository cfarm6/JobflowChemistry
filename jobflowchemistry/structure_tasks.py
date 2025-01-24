from .utils import rdkit2ase, ase2rdkit

from dataclasses import dataclass, field
from jobflow import job, Flow, Maker, Response
from rdkit.Chem import rdDistGeom, rdchem, rdmolops, rdmolfiles
from ase import optimize
from posebusters import PoseBusters
import pickle

@dataclass
class StructureOperation(Maker):
    name: str = "Structure Operation"

    def operation(self, molecule):
        raise NotImplementedError

    @job(files="files", settings="settings", properties="properties")
    def make(self, molecule: rdchem.Mol):
        molecule, properties, settings = self.operation(molecule)
        return Response(
            output={
                "structure": molecule,
                "files": rdmolfiles.MolToV3KMolBlock(molecule),
                "settings": settings,
                "properties": properties,
            },
            stored_data=properties,
        )
