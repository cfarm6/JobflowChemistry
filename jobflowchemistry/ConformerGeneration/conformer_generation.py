from dataclasses import dataclass, field
from typing import Literal, Union
# from dataclasses import dataclass, field
from ..Structure import Structure
from jobflow import job, Flow, Maker, Response

from ..utils import rdkit2ase, ase2rdkit
from rdkit.Chem import rdmolops, rdmolfiles, rdDistGeom, rdchem
import ase.optimize

@dataclass
class ConformerGeneration(Maker):
    name: str = "Conformer Generation"
    def generate_conformers(structure):
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