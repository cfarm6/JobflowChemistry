from dataclasses import dataclass
from typing import Literal, List, Dict

@dataclass
class CRESTLevelBlocks():
    method: Literal[
        "tblite", 
        "gfn2", 
        "gfn1",
        "gfn0",
        "gfnff",
        "orca",
    ] = "gfn2"
    binary: str = "xtb"
    calcspace: str = None
    chrg: int = None
    uhf: int = None
    rdwbo: bool = False
    rddip: bool = False
    dipgrad: bool = False
    gradfile: str = None
    gradtype: Literal["engrad"] = None
    def to_dict(self):
        d = {}
        keys = [
            "method",
            "binary",
            "calcspace",
            "chrg",
            "uhf",
            "rdwbo",
            "rddip",
            "dipgrad",
        ]
        for key, value in vars(self).items():
            if value is None: continue
            d[key] = value
        return d

@dataclass 
class CRESTCalculationBlock():
    name: str = "CREST Calculation Block"
    type: int = 0
    elog: str = None
    eprint: bool = False
    opt_engine: Literal['ancopt', 'rfo', 'gf'] = None
    hess_update: Literal["bfgs", "powell", "sr1", "bofill", "schlegel"] = None
    maxcycle: int = None
    optlev: Literal["crude", "vloose", "loose", "normal", "tight", "vtight", "extreme"] = "normal"
    converge_e: float = None
    converge_g: float = None
    freeze: str = None
    levels: List[CRESTLevelBlocks] = None
    def to_dict(self):
        d = {"type": self.type, "eprint": self.eprint, "optlev": self.optlev}
        if self.elog is not None:
            d['elog'] = self.elog
        if self.opt_engine is not None: 
            d['opt_engine'] = self.opt_engine
        if self.hess_update is not None:
            d['hess_update'] = self.hess_update
        if self.maxcycle is not None:
            d['maxcycle'] = self.maxcycle
        if self.converge_e is not None:
            d['converge_e'] = self.converge_e
        if self.converge_g is not None:
            d['converge_g'] = self.converge_g
        if self.freeze is not None:
            d["freeze"] = self.freeze
        if self.levels is not None:
            d['level'] = list(map(lambda x: x.to_dict(), self.levels))
        return d

@dataclass
class CRESTCalculator():
    name: str = "CREST Calculator"
    threads: int = 1
    runtype: Literal[
        "none",
        "singlepoint",
        "numgrad",
        "optimize",
        "ohess",
        "numhess",
        "ancopt_ensemble",
        "optimize_ensemble",
        "screen_ensemble",
        "ensemble_singlepoints",
        "metadynamics",
        "dynamics",
        "imtd-gc",
        "nci-mtd",
        "entropy",
        "protonate"
    ] = "none"
    blocks: List[Dict] = None
    inputfile: str = None
    def make_dict(self):
        d = {
            'runtype': self.runtype,
            'threads': self.threads,
            'input' : self.inputfile
        }
        for block in self.blocks:
            d[block] = self.blocks[block]
        return d
