from dataclasses import dataclass, field
from typing import Literal, Union
import tomli_w 
import subprocess
import os

from jobflow import job, Flow, Maker, Response

import ase.optimize
from rdkit.Chem import rdmolops, rdmolfiles, rdDistGeom, rdchem

from ..Structure import Structure
from ..utils import rdkit2ase, ase2rdkit

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
                    "settings": {},
                    "properties": [x.output["properties"] for x in jobs],
                },
                addition=jobs,
            )
        if structure.GetNumConformers() > 1:
            jobs = []
            for confId in range(structure.GetNumConformers()):
                s = Structure(rdchem.Mol(structure, confId=confId))
                jobs.append(self.make(s))
            return Response(
                output={
                    "structure": [x.output["structure"] for x in jobs],
                    "settings": {},
                    "properties": [x.output["properties"] for x in jobs],
                },
                addition=jobs,
            )
        structure = self.generate_conformers(structure)
        settings = self.settings()
        properties = self.properties(structure)
        with rdmolfiles.SDWriter("conf.sdf") as f:
            for i in range(structure.GetNumConformers()):
                f.write(structure, confId=i)
        with open("conf.sdf", "r") as f:
            filelines = f.read()
        return Response(output={
            'structure': Structure(structure),
            'files': filelines,
            'settings': settings,
            'properties': properties
        })

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
            if key == "name" or key == "numConfs" or key =="method" or value is None: continue
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
        opts = options(path="input.sdf", **{k: v for k,v in vars(self).items() if v is not None and k != "name"})
        out = main(opts)
        structure.RemoveAllConformers()
        for i, mol in enumerate(rdmolfiles.SDMolSupplier(out, removeHs=False, sanitize=False)):
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
    topo: bool = True
    threads: int = 1
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
    energy_method: Literal[
        "tblite", 
        "gfn2", 
        "gfn1",
        "gfn0",
        "gfnff",
        "orca",
    ] = "gfn2"
    energy_binary: str = "xtb"
    energy_calcspace: str = None
    energy_chrg: int = None
    energy_uhf: int = None
    energy_rdwbo: bool = False
    energy_rddip: bool = False
    energy_dipgrad: bool = False
    energy_gradfile: str = None
    energy_gradtype: Literal["engrad"] = None

    # Metadynamics Block
    dynamics_method: Literal[
        "tblite",
        "gfn2",
        "gfn1",
        "gfn0",
        "gfnff",
        "orca",
    ] = "gfn2"
    dynamics_binary: str = "xtb"
    dynamics_calcspace: str = None
    dynamics_chrg: int = None
    dynamics_uhf: int = None
    dynamics_rdwbo: bool = False
    dynamics_rddip: bool = False
    dynamics_dipgrad: bool = False
    dynamics_gradfile: str = None
    dynamics_gradtype: Literal["engrad"] = None
    def generate_conformers(self, structure: Structure):
        # Write structures to sdf file
        with rdmolfiles.SDWriter("input.sdf") as f:
            f.write(structure)
        self.input = "input.sdf"
        # Empty strucutres
        d = {'calculation':{}, 'cregen': {}}
        calculation_blocks = {"energy":{}, "dynamics": {}}
        # Fill in dictionaries for toml
        for k,v in vars(self).items():
            if v is None: continue
            if k == "name": continue

            if k.split("_")[0] == "energy" or k.split("_")[0] == "dynamics":
                calculation_blocks[k.split("_")[0]][k.split("_")[1]] = v
            elif k in ["ewin", "ethr", "rthr", "bthr"]:
                d['cregen'][k] = v
            else: 
                if k in ["type", "elog", "eprint", "opt_engine","hess_update", "maxcycle","optlev", "converge_e", "converge_g", "freeze"]:
                    d["calculation"][k] = v
                else:
                    d[k] = v
        d["calculation"]["level"] = [calculation_blocks["energy"], calculation_blocks["dynamics"]]
        # Check for charges
        if self.energy_chrg is None or self.dynamics_chrg is None:
            self.energy_chrg = self.dynamics_chrg = rdmolops.GetFormalCharge(structure)
        # Set calculator for dynamics
        d["dynamics"] = {"active": [2]}
        with open("crest.toml", "wb") as f:
            tomli_w.dump(d, f)
        subprocess.call(f"crest --input crest.toml > log.out", shell = True)
        if not os.path.exists("crest_conformers.sdf"): return None
        suppl = rdmolfiles.SDMolSupplier("crest_conformers.sdf", removeHs=False, sanitize=False)
        structure.RemoveAllConformers()
        for m in suppl:
            structure.AddConformer(m.GetConformer(), assignId=True)
        return structure
    def settings(self):
        return {}
    def properties(self, structure):
        return {}
