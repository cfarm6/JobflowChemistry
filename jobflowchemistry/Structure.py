from rdkit.Chem.rdchem import Mol
import pickle

class Structure(Mol):
    periodic: bool = False
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
    """
    Represents a three-dimensional molecular structure.
    Inherits from Structure.
    """
    None

class TwoDimensionalStructure(Structure):
    """
    Represents a two-dimensional molecular structure.
    Inherits from Structure.
    """
    None
