# from pydantic.dataclasses import dataclass
from dataclasses import dataclass, field

# from pydantic import Field
from typing import Literal

from jobflow import job, Flow, Maker, Response
from rdkit.Chem import rdDistGeom, rdmolfiles, rdchem, rdmolops
import subprocess
from posebusters import PoseBusters
from ..Structure import Structure
from ..outputs import Settings, Properties


@dataclass
class StructureGeneration(Maker):
    name: str = field(default="Structure Generation")
    checkStructure: bool = field(default=True)

    def generate_structure(self, molecule: Structure):
        raise NotImplementedError

    @job(files="files", settings="settings", properties="properties")
    def make(self, structure: Structure):
        structure, properties, settings = self.generate_structure(structure)

        buster = PoseBusters(config="mol")
        df = buster.bust(structure, None, None)
        if not (all(df.to_numpy()[0])):
            return Response(stop_children=True)

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
class RDKitGeneration(StructureGeneration):
    name: str = field(default="RDKit Structure Generation")
    method: Literal["ETKDGv3", "ETKDGv2", "ETKDGv1", "ETDG"] = field(
        default="ETKDGv3",
    )

    def generate_structure(self, structure: Structure):

        method = getattr(rdDistGeom, self.method)()
        rdDistGeom.EmbedMolecule(structure, method)
        properties = {}
        settings = {"Structure Generation Method": self.method}
        return structure, properties, settings


@dataclass
class aISSDocking(StructureGeneration):
    name: str = "aISS Docking"
    pocket: bool = False
    no_stack: bool = False
    no_angular: bool = False
    fast: bool = False
    ATM: bool = False
    stepr: float = 2.5
    stepa: float = 45
    maxgen: int = 10
    maxparent: int = 100
    nstack: int = 1000
    nfinal: int = 15
    ensemble: bool = False
    etemp: float = 300.0
    iterations: int = 250
    acc: float = 1.0
    opt: Literal[
        "crude",
        "sloppy",
        "loose",
        "lax",
        "normal",
        "tight",
        "vtight",
        "extreme",
    ] = "normal"
    cycles: int = 100
    optlevel: Literal["gfn1", "gfn2", "gfnff"] = "gfn2"
    alpb: Literal[
        "cetone",
        "acetonitrile",
        "aniline",
        "benzaldehyde",
        "benzene",
        "ch2cl2",
        "chcl3",
        "cs2",
        "dioxane",
        "dmf",
        "dmso",
        "ether",
        "ethylacetate",
        "furane",
        "hexadecane",
        "hexane",
        "methanol",
        "nitromethane",
        "octanol",
        "woctanol",
        "phenol",
        "toluene",
        "thf",
        "water",
    ] = None
    gbsa: Literal[
        "acetone",
        "acetonitrile",
        "benzene",
        "CH2Cl2",
        "CHCl3",
        "CS2",
        "DMF",
        "DMSO",
        "ether",
        "H2O",
        "methanol",
        "n-hexane",
        "THF",
        "toluene",
    ] = None

    @job(files="files", settings="settings", properties="properties")
    def make(self, structure_1: Structure, structure_2: Structure):
        rdmolfiles.MolToMolFile(structure_1, "structure1.mol")
        rdmolfiles.MolToMolFile(structure_2, "structure2.mol")
        charge_1 = rdmolops.GetFormalCharge(structure_1)
        charge_2 = rdmolops.GetFormalCharge(structure_2)
        with open(".CHRG", "w") as f:
            f.write(str(charge_1 + charge_2)+"\n")
            f.write(str(charge_1)+"\n")
            f.write(str(charge_2)+"\n")
        commands = ["xtb", "dock", "structure1.mol", "structure2.mol"]
        if self.pocket:
            commands.append("--pocket")
        if self.no_stack:
            commands.append("--nostack")
        if self.no_angular:
            commands.append("--noangular")
        if self.fast:
            commands.append("--fast")
        if self.ATM:
            commands.append("--atm")
        commands.append(f"--stepr {self.stepr}")
        commands.append(f"--stepa {self.stepa}")
        commands.append(f"--maxgen {self.maxgen}")
        commands.append(f"--maxparent {self.maxparent}")
        commands.append(f"--nstack {self.nstack}")
        commands.append(f"--nfinal {self.nfinal}")
        commands.append(f"--ensemble")
        commands.append(f"--etemp {self.etemp}")
        commands.append(f"--iterations {self.iterations}")
        commands.append(f"--acc {self.acc}")
        commands.append(f"--opt {self.opt}")
        subprocess.call(" ".join(commands), shell=True)
        return 
        # return Response(
        #     output={
        #         "structure": Structure(structure),
        #         "files": rdmolfiles.MolToV3KMolBlock(structure),
        #         "settings": Settings(settings),
        #         "properties": Properties(properties),
        #     },
        #     stored_data=properties,
        # )

