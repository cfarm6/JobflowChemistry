# from pydantic.dataclasses import dataclass
from dataclasses import dataclass, field

# from pydantic import Field
from typing import Literal
import tomli_w
import os
import ase
from icecream import ic

from jobflow import job, Flow, Maker, Response
from rdkit.Chem import rdDistGeom, rdmolfiles, rdchem, rdmolops, rdDetermineBonds
import subprocess
from posebusters import PoseBusters
from ..Structure import Structure
from ..outputs import Settings, Properties
from ..Calculators import CRESTCalculator
from ..utils import ase2rdkit
from ..ConformerGeneration import rdKitConformers


@dataclass
class StructureGeneration(Maker):
    name: str = field(default="Structure Generation")
    check_structure: bool = False

    def generate_structure(self, structure: Structure):
        raise NotImplementedError

    @job(files="files", settings="settings", properties="properties")
    def make(self, structure: Structure):
        if type(structure) is list:
            jobs = [self.make(s) for s in structure]
            return Response(
                output={
                    "structure": [x.output["structure"] for x in jobs],
                    "settings": Settings({}),
                    "files": [x.output["files"] for x in jobs],
                    "properties": [x.output["properties"] for x in jobs],
                },
                detour=jobs,
            )
        if structure.GetNumConformers() > 1:
            jobs = []
            for confId in range(structure.GetNumConformers()):
                s = Structure(rdchem.Mol(structure, confId=confId))
                jobs.append(self.make(s))
            return Response(
                output={
                    "structure": [x.output["structure"] for x in jobs],
                    "files": [x.output["files"] for x in jobs],
                    "settings": Settings({}),
                    "properties": [x.output["properties"] for x in jobs],
                },
                detour=jobs,
            )

        structure = self.generate_structure(structure)
        if type(structure) is list: ic(len(structure))
        if structure is None:
            return Response(stop_children=True)
        settings = self.settings()
        properties = self.properties(structure)

        if self.check_structure:
            buster = PoseBusters(config="mol")
            df = buster.bust(structure, None, None)
            if not (all(df.to_numpy()[0])):
                return Response(stop_children=True)
        if type(structure) is Structure:
            resp = Response(
                output={
                    "structure": Structure(structure),
                    "files": rdmolfiles.MolToV3KMolBlock(structure),
                    "settings": Settings(settings),
                    "properties": Properties(properties),
                },
                stored_data={},
            )
        else:
            resp = Response(
                output={
                    "structure": list(map(lambda x: Structure(x), structure)),
                    "files": list(
                        map(lambda x: rdmolfiles.MolToV3KMolBlock(x), structure)
                    ),
                    "settings": Settings(settings),
                    "properties": Properties(properties),
                },
                stored_data={},
            )
        return resp


@dataclass
class CRESTProtonation(StructureGeneration):
    name: str = "CREST Protonation"
    runtype: Literal["protonate"] = "protonate"
    ion: str = None
    ion_charge: int = 1
    ewin: float = None
    ffopt: bool = True
    freezeopt = None
    finalopt: bool = True
    threads: int = 1

    def make_dict(self):
        keys = ["ion", "ewin", "ffopt", "freezeopt", "finalopt"]
        d = {}
        for k, v in vars(self).items():
            if v is None or k not in keys:
                continue
            d[k] = v
        return d

    def settings(self):
        return {
            "ion": self.ion,
            "ion_charge": self.ion_charge,
            "ewin": self.ewin,
            "ffopt": self.ffopt,
            "freezeopt": self.freezeopt,
            "finalopt": self.finalopt,
            "threads": self.threads,
        }

    def properties(self, structure):
        return {}

    def generate_structure(self, structure: Structure):
        # Write the structure file
        rdmolfiles.MolToXYZFile(structure, "input.xyz")
        # Write the input file
        self.inputfile = "input.xyz"
        d = {"threads": self.threads, "runtype": self.runtype, "input": "input.xyz"}
        d["protonation"] = self.make_dict()
        with open("crest.toml", "wb") as f:
            tomli_w.dump(d, f)
        # Run the calculation
        chrg = rdmolops.GetFormalCharge(structure)
        subprocess.call(f"crest --input crest.toml {chrg} > log.out", shell=True)
        # Parse the output
        if not os.path.exists("protonated.xyz"):
            return Response(stop_children=True)
        # Generate bonds with obabel
        subprocess.call("obabel protonated.xyz -O protonated.sdf", shell=True)
        # Return array of adducts
        suppl = list(
            rdmolfiles.SDMolSupplier("protonated.sdf", sanitize=False, removeHs=False)
        )
        if len(suppl) == 1:
            suppl = suppl[0]

        return suppl


@dataclass
class CRESTDeprotonation(
    StructureGeneration,
):
    name: str = "CREST Deprotonation"
    ewin: float = None

    def generate_structure(self, structure: Structure):
        with rdmolfiles.SDWriter("input.sdf") as f:
            f.write(structure)
        d = {"calculation": {"level": {"method": "gfn2", "rdwbo": True}}}
        with open("crest.toml", "wb") as f:
            tomli_w.dump(d, f)
        commands = ["crest", "input.sdf", "--deprotonate", "--input", "crest.toml"]
        if self.ewin is not None:
            commands.append(f"--ewin {self.ewin}")
        charge = rdmolops.GetFormalCharge(structure)
        commands.append(f"--chrg {charge} --newversion")
        commands.append(" > log.out")
        subprocess.call(" ".join(commands), shell=True)
        if not os.path.exists("deprotonated.xyz"):
            return None
        subprocess.call("obabel deprotonated.xyz -O deprotonated.sdf", shell=True)

        deprotonated_conformers = list(
            rdmolfiles.SDMolSupplier(
                "deprotonated.sdf", 
                sanitize=False, 
                removeHs=False, 
                strictParsing=False
            )
        )

        with rdmolfiles.SDWriter(
            "updated.sdf",
        ) as f:
            for i, s in enumerate(deprotonated_conformers):
                rdDetermineBonds.DetermineBonds(s, charge=-1)
                # topo = rdchem.Mol(structure) # Structure with only single bonds
                # _s = rdchem.Mol(s)
                # ic(i)
                # for mol in [_s, topo]:
                #     for bond in mol.GetBonds():
                #         bond.SetStereo(rdchem.BondStereo.STEREOANY)
                #         bond.SetBondType(rdchem.BondType.SINGLE)
                # idxs = topo.GetSubstructMatch(_s)  # s_idx -> structure_idx
                # # Find the missing atom idx
                # for atom in topo.GetAtoms():
                #     if atom.GetIdx() in idxs:
                #         continue
                #     h_idx = atom.GetIdx()
                #     bond = atom.GetBonds()[0]
                #     begin_idx = bond.GetBeginAtomIdx()
                #     end_idx = bond.GetEndAtomIdx()
                #     other_idx = begin_idx
                #     if other_idx == h_idx:
                #         other_idx = end_idx
                #     ic(other_idx)
                #     ic(idxs)
                #     s_idx = idxs.index(other_idx)
                #     s_atom = s.GetAtomWithIdx(s_idx)
                #     s_atom.SetFormalCharge(-1)
                f.write(s)
        structures = list(rdmolfiles.SDMolSupplier("updated.sdf", removeHs=False, sanitize=False))
        if len(structures) == 1:
            return Structure(structures[0])
        return structures

    def settings(self):
        return {"ewin": self.ewin}

    def properties(self, structure):
        return {}


@dataclass
class RDKitGeneration(StructureGeneration):
    name: str = "rdKit Generation"
    method: Literal["ETDG", "ETKDG", "ETKDGv2", "ETKDGv3", "KDG", "srETKDGv3"] = (
        "ETKDGv3"
    )
    boundsMatForceScaling: float = None
    boxSizeMult: float = None
    clearConfs: bool = None
    embedFragmentsSeparately: bool = None
    enableSequentialRandomSeeds: bool = None
    enforceChirality: bool = None
    forceTransAmides: bool = None
    ignoreSmoothingFailures: bool = None
    maxIterations: int = None
    numThreads: int = 1
    numZeroFail: int = None
    onlyHeavyAtomsForRMS: bool = None
    optimizerForceTol: float = None
    pruneRmsThresh: float = None
    randNegEig: bool = None
    randomSeed: int = None
    symmetrizeConjugatedTerminalGroupsForPruning: bool = None
    trackFailures: bool = None
    useBasicKnowledge: bool = None
    useExpTorsionAnglePrefs: bool = None
    useMacrocycle14config: bool = None
    useMacrocycleTorsions: bool = None
    useRandomCoords: bool = None
    useSmallRingTorsions: bool = None
    useSymmetryForPruning: bool = None

    def generate_structure(self, structure: Structure):
        params = getattr(rdDistGeom, self.method)()
        for key, value in vars(self).items():
            if key == "name" or key == "method" or value is None:
                continue
            setattr(params, key, value)
        rdDistGeom.EmbedMultipleConfs(structure, 1, params)
        return structure

    def settings(self):
        d = {}
        for key, value in vars(self).items():
            if key == "name" or value is None:
                continue
            d[key] = value
        return d

    def properties(self, structure: Structure):
        return {}


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
        if type(structure_1) is list or type(structure_2) is list:
            if type(structure_1) is not list:
                structure_1 = [structure_1]
            if type(structure_2) is not list:
                structure_2 = [structure_2]
            jobs = []
            for s1 in structure_1:
                for s2 in structure_2:
                    jobs.append(self.make(s1, s2))
            return Response(replace=jobs)
        rdmolfiles.MolToMolFile(structure_1, "structure1.mol")
        rdmolfiles.MolToMolFile(structure_2, "structure2.mol")
        charge_1 = rdmolops.GetFormalCharge(structure_1)
        charge_2 = rdmolops.GetFormalCharge(structure_2)
        with open(".CHRG", "w") as f:
            f.write(str(charge_1 + charge_2) + "\n")
            f.write(str(charge_1) + "\n")
            f.write(str(charge_2) + "\n")
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
        if self.ensemble:
            commands.append(f"--ensemble")
        commands.append(f"--etemp {self.etemp}")
        commands.append(f"--iterations {self.iterations}")
        commands.append(f"--acc {self.acc}")
        commands.append(f"--opt {self.opt}")
        commands.append(" > log.out")
        subprocess.call(" ".join(commands), shell=True)
        with open("run.sh", "w") as f:
            f.write(" ".join(commands))
        if not os.path.exists("final_structures.xyz"):
            return Response(stop_children=True)
        subprocess.call(
            "obabel final_structures.xyz -O final_structures.sdf", shell=True
        )
        suppl = list(
            rdmolfiles.SDMolSupplier(
                "final_structures.sdf", sanitize=False, removeHs=False
            )
        )
        if len(suppl) == 1:
            resp = Response(
                output={
                    "structure": Structure(suppl[0]),
                    "files": rdmolfiles.MolToV3KMolBlock(suppl[0]),
                    "settings": Settings({}),
                    "properties": Properties({}),
                },
                stored_data={},
            )
        else:
            resp = Response(
                output={
                    "structure": list(map(lambda x: Structure(x), suppl)),
                    "files": list(map(lambda x: rdmolfiles.MolToV3KMolBlock(x), suppl)),
                    "settings": Settings({}),
                    "properties": Properties({}),
                },
                stored_data={},
            )
        return resp
