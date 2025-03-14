# from pydantic.dataclasses import dataclass
from dataclasses import dataclass
from jobflow import job, Response, Maker

from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdchem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import AddHs
from typing import Union
import requests
import pickle
from ..Structure import Structure
from ..outputs import Settings, Properties
from ..inputs import SMARTS, SMILES, PubChemCID

@dataclass
class StructureInput(Maker):
    name: str = "Input to Molecule"
    removeSalts: bool = True

    def get_structure(self, input):
        raise NotImplementedError

    def make_structure(self, input):
        structure, properties, settings = self.get_structure(input)
        if self.removeSalts:
            structure = SaltRemover().StripMol(structure)
        structure = AddHs(structure)
        return Response(
            output={
                "structure": Structure(structure),
                "files": rdmolfiles.MolToV3KMolBlock(structure),
                "settings": Settings(settings),
                "properties": Properties(properties),
            },
            stored_data=properties,
        )
    @job(files="files", settings="settings", properties="properties")
    def make(self, input: Union[int, str]):
        raise NotImplementedError

@dataclass
class PubChemInput(StructureInput):
    name: str = "PubChem CID To Graph"

    def get_structure(self, input: int):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{input}/sdf"
        resp = requests.post(url)
        with open("log", "w") as f:
            f.write(resp.content.decode("utf-8"))
        molecule = rdmolfiles.MolFromMolBlock(resp.content.decode("utf-8"), removeHs=False)
        properties = {}
        settings = {'Input Method': 'PubChem CID'}
        return molecule, properties, settings

    @job(files="files", settings="settings", properties="properties")
    def make(self, input: int):
        return super().make_structure(input)


@dataclass
class SmilesInput(StructureInput):
    name: str = "SMILES To Graph"

    def get_structure(self, smiles: str):
        molecule = rdmolfiles.MolFromSmiles(smiles)
        print(len(list(molecule.GetAtoms())))
        properties = {}
        settings = {'Input Method': 'SMILES Input'}
        return molecule, properties, settings

    @job(files="files", settings="settings", properties="properties")
    def make(self, SMILES: SMILES):
        return super().make_structure(SMILES)


@dataclass
class SmartsInput(StructureInput):
    name: str = "SMARTS To Graph"

    def get_structure(self, smarts: str):
        molecule = rdmolfiles.MolFromSmarts(smarts)
        properties = {}
        settings = {'Input Method': 'SMARTS Input'}
        return molecule, properties, settings

    @job(files="files", settings="settings", properties="properties")
    def make(self, SMARTS: SMARTS):
        return super().make_structure(SMARTS)
