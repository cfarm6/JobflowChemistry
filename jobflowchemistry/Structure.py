from rdkit.Chem.rdchem import Mol
import pickle

class Structure(Mol):
    periodic: bool = False
    def as_dict(self):
        # pickle the molecule instance, not super() — the latter unpickles as
        # bare `super` and breaks RDKit C++ calls (e.g. EmbedMultipleConfs).
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "data": pickle.dumps(self),
        }
        return d

    @classmethod
    def from_dict(cls, d):
        if type(d["data"]) is str:
            obj = pickle.loads(eval(d["data"]))
        else:
            obj = pickle.loads(d["data"])
        if type(obj) is super:
            obj = obj.__self__
        if isinstance(obj, cls):
            return obj
        if isinstance(obj, Mol):
            return cls(obj)
        raise TypeError(f"Cannot restore {cls.__name__} from {type(obj)!r}")

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
