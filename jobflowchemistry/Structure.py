from rdkit.Chem.rdchem import Mol
import pickle
from icecream import ic

class Structure(Mol):
    def as_dict(self):
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "data": pickle.dumps(super()),
        }
        return d
    @classmethod
    def from_dict(cls, d):
        if type(d["data"]) is str:
            return pickle.loads(eval(d["data"]))
        else: 
            return pickle.loads(d["data"])

class ThreeDimensionalStructure(Structure):
    None

class TwoDimensionalStructure(Structure):
    None
