from dataclasses import dataclass, fields
from jobflow import Response, job, Maker
import subprocess

# - - - - - - - -
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from typing import Union, Literal, List
from io import StringIO
# - - - - - - - -
from ..Structure import Structure
from ..outputs import Settings, Properties


@dataclass
class LigandDocking(Maker):
    name: str = "Ligand Docking"

    def dock(self, ligand, protein):
        raise NotImplementedError

    def get_settings(self):
        raise NotImplementedError

    def get_properties(self, docked_structures, energies):
        raise NotImplementedError

    @job(files="files", settings="settings", properties="properties")
    def make(self, structure: Structure, protein):
        if type(structure) is list:
            jobs = [self.make(s, protein) for s in structure]
            return Response(
                output={
                    "structure": [x.output["structure"] for x in jobs],
                    "settings": Settings({}),
                    "properties": [x.output["properties"] for x in jobs],
                },
                replace=jobs,
            )
        if structure.GetNumConformers() > 1:
            jobs = []
            for confId in range(structure.GetNumConformers()):
                s = Structure(rdchem.Mol(structure, confId=confId))
                jobs.append(self.make(s, protein))
            return Response(
                output={
                    "structure": [x.output["structure"] for x in jobs],
                    "settings": Settings({}),
                    "properties": [x.output["properties"] for x in jobs],
                },
                replace=jobs,
            )

        structures, energies = self.dock(structure, protein)
        settings = self.get_settings()
        properties = self.get_properties(structures, energies)
        
        sio = StringIO()
        writer = rdmolfiles.SDWriter(sio)
        cids = list(map(lambda x: x.GetId(), structure.GetConformers()))
        for cid in cids:
            writer.write(structures, confId=cid)
        writer.flush()
        struct_file = sio.getvalue()
        return Response(
            output={
                "structure": Structure(structure),
                "files": struct_file,
                "settings": Settings(settings),
                "properties": Properties(properties),
            },
            stored_data=properties,
        )


@dataclass
class AutoDockVina(LigandDocking):
    name: str = "AutoDock Vina Ligand Docking"
    scoring_function: Literal["vina", "vinardo", "ad4"] = "vina"
    cpu: int = 0
    seed: int = 0
    no_refine: bool = False
    exhaustiveness: int = 8
    n_poses: int = 20
    min_rmsd: float = 1.0
    max_evals: int = 0
    center: List[float] = None
    box_size: List[float] = None
    space: float = 0.375
    energy_range: float = 3.0
    force_even_voxels: bool = False

    def dock(self, ligand: Structure, protein: str):
        from meeko import MoleculePreparation
        from meeko import PDBQTWriterLegacy
        from meeko import PDBQTMolecule
        from meeko import RDKitMolCreate
        from meeko import Polymer

        from vina import vina

        # Prep Ligand
        mk_prep = MoleculePreparation()
        molsetup_list = mk_prep(rdchem.Mol(ligand))
        molsetup = molsetup_list[0]
        ligand_pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)[0]
        print(ligand_pdbqt_string)
        # Prep Protein
        with open("input.pdb", "w") as f:
            f.write(protein)
        # protein = Polymer.from_pdb_string(protein, chem_templates=None, mk_prep=T)
        # pdbqt_string, _ = PDBQTWriterLegacy.write_string_from_polymer(protein)
        subprocess.run(
            "mk_prepare_receptor.py -i input.pdb -o prepared -p ", shell=True
        )
        # with open("prepared.pdbqt", "r") as f:
        #     protein_pdbqt_string = f.read()

        # Initialize Docking
        vina = vina.Vina(
            sf_name=self.scoring_function,
            cpu=self.cpu,
            seed=self.seed,
            no_refine=self.no_refine,
        )
        vina.set_receptor("prepared.pdbqt")
        vina.set_ligand_from_string(ligand_pdbqt_string)
        vina.compute_vina_maps(
            center=self.center,
            box_size=self.box_size,
            spacing=self.space,
            force_even_voxels=self.force_even_voxels,
        )
        vina.dock(
            exhaustiveness=self.exhaustiveness,
            n_poses=self.n_poses,
            min_rmsd=self.min_rmsd,
            max_evals=self.max_evals,
        )
        poses = vina.poses(n_poses=self.n_poses, energy_range=self.energy_range)
        pdbqt_molecule = PDBQTMolecule(
            poses,
            skip_typing=True,
        )
        rdmolList = RDKitMolCreate.from_pdbqt_mol(
            pdbqt_molecule,
        )[0]
        energies = vina.energies(n_poses=self.n_poses, energy_range=self.energy_range)
        return Structure(rdmolList), energies

    def get_properties(self, docked_structures, energies):
        fields = ["Total", "Inter", "Intra", "Torsions", "Intra best pose"]
        if self.scoring_function == "ad4":
            fields = ["Total", "Inter", "Intra", "Torsions", "-Intra"]
        results = []
        for i, energy in enumerate(energies):
            results.append({x[0]: x[1] for x in zip(fields, energy)})
        return {"Global": {f"{self.scoring_function} Docking Scores": results}}

    def get_settings(self):
        return {
            "scoring function": self.scoring_function,
            "cpu": self.cpu,
            "seed": self.seed,
            "no_refine": self.no_refine,
            "exhaustiveness": self.exhaustiveness,
            "n_poses": self.n_poses,
            "min_rmsd": self.min_rmsd,
            "max_evals": self.max_evals,
            "center": self.center,
            "box_size": self.box_size,
            "space": self.space,
            "energy_range": self.energy_range,
            "force_even_voxels": self.force_even_voxels,
        }
