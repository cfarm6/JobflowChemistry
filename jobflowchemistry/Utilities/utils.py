from dataclasses import dataclass, field
from jobflow import Response, job, Flow, Maker
import numpy as np
from rdkit.Chem import rdmolfiles
from icecream import ic
from ..Structure import Structure
from ..outputs import Properties
from typing import Literal

@dataclass
class BoltzmannWeighting(Maker):
    name: str = "Boltzmann Weighting"
    temperature: float = 300.0
    property_type: Literal["global"] = "global"
    property: str = None
    energy_name: str = "Total Energy [eV]"

    def get_settings(self):
        return {k:v for k,v in vars(self).items() if k != "name"}
    
    @job(files="files", settings="settings", properties="properties")
    def make(self, structure: Structure):
        # molecule = pickle.loads(molecule)
        if type(structure) is not list: Response(stop_children=True)

        energy = [s.GetProp(self.energy_name, autoConvert=True) for s in structure]
        energy = np.array(energy)
        energy = energy - np.min(energy)

        if self.property_type == "global":
            prop = [s.GetProp(self.property, autoConvert=True) for s in structure]
        
        prop = np.array(prop)
        kb = 8.617333262e-5  # eV/K
        num = np.sum(prop * np.exp(-1 * energy / kb / self.temperature))
        denom = np.sum(np.exp(-1 * energy / kb / self.temperature))
        average = num / denom
        ic(prop)
        ic(average)
        settings = self.get_settings()

        properties = {self.property_type: {f"Boltzmann Weighted {self.property}": average}}
        
        return Response(
            output={
                "structure": [Structure(s) for s in structure],
                "files": [rdmolfiles.MolToV3KMolBlock(s) for s in structure],
                "settings": settings,
                "properties": Properties(properties),
            },
            stored_data=properties,
        )
