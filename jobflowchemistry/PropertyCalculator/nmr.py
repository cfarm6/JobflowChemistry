from .property_calculator import PropertyCalculator

from dataclasses import dataclass, field
from jobflow import Response, job, Flow, Maker

@dataclass
class xTBNMR(PropertyCalculator):
    name: str = "xTB NMR Prediction"