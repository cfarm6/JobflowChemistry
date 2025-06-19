from dataclasses import dataclass, field
from typing import Literal, Union
import tomli_w
import subprocess
import os

from jobflow import job, Flow, Maker, Response

from rdkit.Chem import rdmolops, rdmolfiles, rdDistGeom, rdchem

from ..Structure import Structure
from ..utils import rdkit2ase, ase2rdkit
from ..outputs import Settings


@dataclass
class ConformerGeneration(Maker):
    name: str = "Conformer Generation"

    def generate_conformers(self, structure):
        raise NotImplementedError

    def settings(self):
        raise NotImplementedError

    def properties(self, structure):
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
                replace=jobs,
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
                replace=jobs,
            )
        structure = self.generate_conformers(structure)
        settings = self.settings()
        properties = self.properties(structure)
        with rdmolfiles.SDWriter("conf.sdf") as f:
            for i in range(structure.GetNumConformers()):
                f.write(structure, confId=i)
        with open("conf.sdf", "r") as f:
            filelines = f.read()
        return Response(
            output={
                "structure": Structure(structure),
                "files": filelines,
                "settings": settings,
                "properties": properties,
            }
        )


@dataclass
class rdKitConformers(ConformerGeneration):
    name: str = "rdKit Conformers"
    method: Literal["ETDG", "ETKDG", "ETKDGv2", "ETKDGv3", "KDG", "srETKDGv3"] = (
        "ETKDGv3"
    )
    numConfs: int = 1
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

    def generate_conformers(self, structure: Structure):
        params = getattr(rdDistGeom, self.method)()
        for key, value in vars(self).items():
            if key == "name" or key == "numConfs" or key == "method" or value is None:
                continue
            setattr(params, key, value)
        rdDistGeom.EmbedMultipleConfs(structure, self.numConfs, params)
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
class Auto3DConformers(ConformerGeneration):
    name: str = "Auto 3D Conformer Generation"
    k: int = None
    window: float = None
    enumerate_tautomer: bool = False
    enumerate_isomer: bool = True
    mpi_np: int = 4
    max_confs: int = None
    use_gpu: bool = True
    gpu_idx: int = 0
    capacity: int = 42
    optimizing_engine: Literal["ANI2x", "ANI2xt", "AIMNET"] = "AIMNET"
    patience: int = 1000
    opt_steps: int = 5000
    convergence_threshold: float = 0.003
    threshold: float = 0.3
    memory: int = None
    batchsize_atoms: int = 1024

    def generate_conformers(self, structure: Structure):
        from Auto3D.auto3D import options, main

        structure.SetProp("_Name", "auto3d")
        with rdmolfiles.SDWriter("input.sdf") as f:
            f.write(structure)
        opts = options(
            path="input.sdf",
            **{k: v for k, v in vars(self).items() if v is not None and k != "name"},
        )
        out = main(opts)
        structure.RemoveAllConformers()
        for i, mol in enumerate(
            rdmolfiles.SDMolSupplier(out, removeHs=False, sanitize=False)
        ):
            structure.AddConformer(mol.GetConformer(), assignId=True)
        return structure

    def settings(self):
        return {k: v for k, v in vars(self).items() if v is not None}

    def properties(self, structure):
        return {}


@dataclass
class CRESTConformers(ConformerGeneration):
    name: str = "CREST Conformer Generation"
    runtype: Literal["imtd-gc", "nci-mtd", "imtd-smtd"] = "imtd-gc"
    preopt: bool = True
    multilevelopt: bool = True
    topo: bool = True
    parallel: int = 1
    opt_engine: Literal["ancopt", "rfo", "gd"] = "ancopt"
    hess_update: Literal["bfgs", "powell", "sd1", "bofill", "schlegel"] = "bfgs"
    maxcylcle: int = None
    optlev: Literal[
        "crude", "vloose", "loose", "normal", "tight", "vtight", "extreme"
    ] = "normal"
    converge_e: float = None
    converge_g: float = None
    freeze: str = None
    # CREGEN Block
    ewin: float = 6.0
    ethr: float = 0.05
    rthr: float = 0.125
    bthr: float = 0.01
    # Optimization Calculation Block
    calculation_energy_method: Literal[
        "tblite",
        "gfn2",
        "gfn1",
        "gfn0",
        "gfnff",
        "orca",
    ] = "gfn2"
    calculation_energy_binary: str = "xtb"
    calculation_energy_calcspace: str = None
    calculation_energy_chrg: int = None
    calculation_energy_uhf: int = None
    calculation_energy_rdwbo: bool = False
    calculation_energy_rddip: bool = False
    calculation_energy_dipgrad: bool = False
    calculation_energy_gradfile: str = None
    calculation_energy_gradtype: Literal["engrad"] = None

    # Metadynamics Block
    calculation_dynamics_method: Literal[
        "tblite",
        "gfn2",
        "gfn1",
        "gfn0",
        "gfnff",
        "orca",
    ] = "gfn2"
    calculation_dynamics_binary: str = "xtb"
    calculation_dynamics_calcspace: str = None
    calculation_dynamics_chrg: int = None
    calculation_dynamics_uhf: int = None
    calculation_dynamics_rdwbo: bool = False
    calculation_dynamics_rddip: bool = False
    calculation_dynamics_dipgrad: bool = False
    calculation_dynamics_gradfile: str = None
    calculation_dynamics_gradtype: Literal["engrad"] = None
    dynamics_dump: float = 100.0

    def generate_conformers(self, structure: Structure):

        # Write structures to sdf file
        with rdmolfiles.SDWriter("input.sdf") as f:
            f.write(structure)
        if self.calculation_energy_chrg is None:
            self.calculation_energy_chrg = rdmolops.GetFormalCharge(structure)
        if self.calculation_dynamics_chrg is None:
            self.calculation_dynamics_chrg = rdmolops.GetFormalCharge(structure)
        self.input = "input.sdf"
        # Empty strucutres
        d = {"calculation": {}, "cregen": {}, "dynamics": {"active": [2]}}
        calculation_blocks = {"energy": {}, "dynamics": {}}
        # Fill in dictionaries for toml
        for k, v in vars(self).items():
            if v is None:
                continue
            if k == "name":
                continue
            if k.split("_")[0] == "calculation":
                if k.split("_")[1] == "energy" or k.split("_")[1] == "dynamics":
                    calculation_blocks[k.split("_")[1]][k.split("_")[2]] = v
            elif k.split("_")[0] == "dynamics":
                d["dynamics"][k.split("_")[1]] = v
            elif k in ["ewin", "ethr", "rthr", "bthr"]:
                d["cregen"][k] = v
            else:
                if k in [
                    "type",
                    "elog",
                    "eprint",
                    "opt_engine",
                    "hess_update",
                    "maxcycle",
                    "optlev",
                    "converge_e",
                    "converge_g",
                    "freeze",
                ]:
                    d["calculation"][k] = v
                else:
                    d[k] = v
        d["calculation"]["level"] = [
            calculation_blocks["energy"],
            calculation_blocks["dynamics"],
        ]
        # Check for charges
        if (
            self.calculation_energy_chrg is None
            or self.calculation_dynamics_chrg is None
        ):
            self.calculation_energy_chrg = self.calculation_dynamics_chrg = (
                rdmolops.GetFormalCharge(structure)
            )

        with open("crest.toml", "wb") as f:
            tomli_w.dump(d, f)
        subprocess.call(f"crest --input crest.toml --noreftopo --mquick > log.out", shell=True)
        if not os.path.exists("crest_conformers.sdf"):
            return None
        suppl = rdmolfiles.SDMolSupplier(
            "crest_conformers.sdf", removeHs=False, sanitize=False
        )
        structure.RemoveAllConformers()
        for m in suppl:
            structure.AddConformer(m.GetConformer(), assignId=True)
        return structure

    def settings(self):
        return {}

    def properties(self, structure):
        return {}
