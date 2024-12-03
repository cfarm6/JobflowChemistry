from dataclasses import dataclass, field
from typing import Literal
# from dataclasses import dataclass, field
from ..structure_tasks import StructureOperation
from ..utils import rdkit2ase, ase2rdkit
from rdkit.Chem import rdmolops
import ase.optimize

@dataclass
class ConformerGeneration(StructureOperation):
    name: str = "Conformer Generation"

    method: Literal["CREST", "Auto3d"] = field(default="Auto3d")

    window: float = 10.0  # kcal/mol

    def process_out(self, output):
        raise NotImplementedError

    def operation(self, molecule):
        properties = {}
        if self.method not in self.methods:
            raise ValueError

        return properties, molecule
