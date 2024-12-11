# from pydantic import Field
from dataclasses import dataclass, field

# from pydantic.dataclasses import dataclass

from jobflow import Response, job, Flow, Maker
import pickle
from rdkit.Chem import rdmolops, rdmolfiles, rdchem
from rdkit import Chem
from typing import Any, Tuple, Union, Dict
import os
from icecream import ic
from typing import Literal
import subprocess

from ..EnergyCalculation import EnergyCalculation
from ..Structure import Structure
from ..outputs import Properties, Settings


@dataclass
class PropertyCalculator(Maker):
    name: str = "Property Calculator"

    def get_settings(self):
        raise NotImplementedError

    def get_properties(self, molecule: Structure):
        raise NotImplementedError

    @job(files="files", settings="settings", properties="properties")
    def make(self, structure: Structure):
        # molecule = pickle.loads(molecule)
        if type(structure) is list:
            ic(len(structure))
            jobs = [self.make(s) for s in structure]
            return Response(
                output={
                    "structure": [x.output["structure"] for x in jobs],
                    "settings": Settings({}),
                    "files": [x.output["files"] for x in jobs],
                    "properties": [x.output["properties"] for x in jobs],
                },
                addition=jobs,
            )
        if structure.GetNumConformers() > 1:
            jobs = []
            ic(structure.GetNumConformers())
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
                addition=jobs,
            )
        resp = self.get_properties(structure)
        ic(resp)
        structure = resp[0]
        properties = resp[1]
        if "global" in properties:
            for k,v in properties["global"].items():
                structure.SetDoubleProp(k, float(v), computed=True)
        if "atomic" in properties:
            for k,v in properties["atomic"].items():
                for i, atom in enumerate(structure.GetAtoms()):
                    atom.SetDoubleProp(k, float(v[i]))
        if "bond" in properties:
            for k,v in properties["bond"].items():
                bond = structure.GetBondBetweenAtoms(v[0], v[1])
                bond.SetDoubleProp(k, float(v[2]))

        settings = self.get_settings()
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
class CollisionCrossSectionCalculator(PropertyCalculator):
    name: str = "Collision Cross Section Calculator"
    pass


@dataclass
class MobcalCollisionCrossSection(CollisionCrossSectionCalculator):
    name: str = "Mobcal CCS Calculation"
    charge_type: str = None
    I2: int = 123456
    buffer_gas: Literal["nitrogen"] = "nitrogen"
    buffer_gas_mass: float = 28.014
    temperature: float = 300.0
    ipr: int = 1000
    itn: int = 10
    inp: int = 40
    imp: int = 25
    num_threads: int = 1
    def get_settings(self):
        return {k:v for k,v in vars(self).items() if k != "name"}
    def get_properties(self, structure: Structure):
        if self.charge_type not in structure.GetAtomWithIdx(0).GetPropNames():
            ic(self.charge_type)
            ic(list(structure.GetPropNames()))
            return Response(stop_children=True)
        with open("atomtype_parameters.in", "w") as f:
            f.write(
                """ATOM    MASS        EOLJ        ROLJ    RHS
                    H   1.0079      0.0189      1.2409  2.20
                    C   12.0107     0.0977      3.5814  2.70
                    N   14.0067     0.0828      4.3920  2.70
                    O   15.9994     0.0558      3.2550  2.70
                    F   18.9984     0.0465      3.1285  2.70
                    Na+ 22.9897     0.0300      2.9830  2.853
                    Si  28.0855     0.4020      4.2950  2.95
                    P   30.9738     0.305       4.1470  4.20
                    S   32.065      0.2740      4.0350  3.50
                    Cl  35.453      0.0465      3.1285  2.70
                    Fe  55.845      0.0130      2.9120  3.50"""
            )
        with open("mobcal.params", "w") as f:
            for k,v in vars(self).items():
                if k not in ["I2", "buffer_gas", "buffer_gas_mass", "temperature", "ipr", "itn", "inp", "imp" , "num_threads"]: continue
                f.write(f"{k.upper()} {v}\n")
        positions = structure.GetConformer().GetPositions()
        masses = [atom.GetMass() for atom in structure.GetAtoms()]
        charges = [atom.GetProp(self.charge_type) for atom in structure.GetAtoms()]

        with open("mobcal.mfj", "w") as f:
            f.write(self.name + "\n")
            f.write("1\n")
            f.write(str(structure.GetNumAtoms())+"\n")
            f.write("ang\n")
            f.write("calc\n")
            f.write("1.0")
            for position, mass, charge in zip(positions, masses, charges):
                f.write(f"\n{position[0]} {position[1]} {position[2]} {mass} {charge}")

        subprocess.call(
            "mobcal_shm mobcal.params atomtype_parameters.in mobcal.mfj log.out",
            shell=True,
        )
        CCS = 0.0
        with open("log.out", "r") as f:
            lines = f.readlines()
            mobility_line = list(
                filter(
                    lambda x: "inverse average (second order) TM mobility" in x, lines
                )
            )[-1]
            mobility = float(mobility_line.split()[-1])

            CCS_line = list(filter(lambda x: "average TM cross section" in x, lines))[
                -1
            ]
            CCS = float(CCS_line.split()[-1])

            standard_deviation_percent_line = list(
                filter(lambda x: "standard deviation (percent)" in x, lines)
            )[-1]
            standard_deviation_percent = float(
                standard_deviation_percent_line.split()[-1]
            )
            standard_deviation = standard_deviation_percent / 100 * CCS
        properties = {
            "global": {
                "CCS [A^2]": CCS,
                "CCS Standard Deviation [A^2]": standard_deviation,
                "Mobility": 1 / mobility,
            }
        }
        return structure, properties


@dataclass
class QupKakePrediction(PropertyCalculator):
    name: str = "QupKake pKa Prediction"
    mp: Union[bool, int] = True
    tautomerize: bool = False

    def get_settings(self):
        return {"mp": self.mp, "tautomerize": self.tautomerize}

    def get_properties(self, structure: Structure):
        from qupkake import predict as qpka
        from qupkake.mol_dataset import MolDataset

        if not os.path.exists("raw"):
            os.mkdir("raw")
        if not os.path.exists("processed"):
            os.mkdir("processed")
        structure.SetProp("_Name", self.name)
        with Chem.SDWriter("raw/mol.sdf") as f:
            f.write(structure)
        filename = "mol.sdf"
        mol_col = "ROMol"
        name_col = "name"
        root = "."
        output = "output.sdf"
        dataset = MolDataset(
            root=root,
            filename=filename,
            tautomerize=self.tautomerize,
            name_col=name_col,
            mol_col=mol_col,
            mp=self.mp,
        )

        prot_model, deprot_model, pka_model = qpka.load_models()
        prot_indices = qpka.predict_sites(dataset, prot_model)
        deprot_indices = qpka.predict_sites(dataset, deprot_model)
        qpka.make_sites_prediction_files(
            root, dataset, prot_indices, deprot_indices, output
        )
        pair_dataset = qpka.load_mol_pair_dataset(
            root=root,
            filename=output,
            name_col=name_col,
            mol_col="ROMol",
            idx_col="idx",
            type_col="pka_type",
            mp=self.mp,
        )
        pka_predictions = qpka.predict_pka(pair_dataset, pka_model)
        return structure, {"global": {"pKa": pka_predictions.tolist()}}


@dataclass
class CalculateGlobalConceptualDFTProperties(PropertyCalculator):
    name: str = "Calculate Global CDFT Properties"

    @job(files="files", settings="settings", properties="properties")
    def make(
        self,
        properties_neutral: Properties,
        properties_cation: Properties,
        properties_anion: Properties,
    ):
        neutral_state_energy = properties_neutral["global"]["Total Energy [eV]"]
        cation_state_energy = properties_cation["global"]["Total Energy [eV]"]
        anion_state_energy = properties_anion["global"]["Total Energy [eV]"]
        vertical_ionization_potential = anion_state_energy - neutral_state_energy
        vertical_eletron_affinity = neutral_state_energy - cation_state_energy
        mulliken_electronegativity = (
            vertical_ionization_potential + vertical_eletron_affinity
        ) / 2
        chemical_potential = -mulliken_electronegativity
        hardness = vertical_ionization_potential - vertical_eletron_affinity
        softness = 1 / hardness
        elctrophilicity_index = chemical_potential**2 / 2 / hardness

        properties = {
            "global": {
                "Vertical Ionization Potential [eV]": vertical_ionization_potential,
                "Vertical Eletron Affinity [eV]": vertical_eletron_affinity,
                "Mulliken Electronegativity [eV]": mulliken_electronegativity,
                "Chemical Potential [eV]": chemical_potential,
                "Hardness [eV]": hardness,
                "Softness [1/eV]": softness,
                "Electrophilicity Index [eV]": elctrophilicity_index,
            }
        }
        return Response(
            output={
                "properties": Properties(properties),
            },
            stored_data=properties,
        )
