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

        structure, properties = self.get_properties(structure)

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


@dataclass
class ConceptualDFTWorkflow(PropertyCalculator):
    name: str = "Conceptual DFT Workflow"

    @job(files="files", settings="settings", properties="properties")
    def make(
        self,
        calculator: EnergyCalculation,
        structure: Structure,
        properties_neutral: Properties = None,
    ):
        if type(structure) is list:
            jobs = [self.make(s) for s in structure]
            return Response(replace=jobs)
        jobs = []
        if calculator.charge is None:
            calculator.charge = rdmolops.GetFormalCharge(structure)
        baseCharge = calculator.charge
        if properties_neutral is None:
            calculator.charge = baseCharge
            N = calculator.make(structure=structure)
            jobs.append(N)
            properties_neutral = N.output["properties"]
        N_energy = properties_neutral

        calculator.charge = baseCharge + 1
        Np1 = calculator.make(structure=structure)
        Np1_energy = Np1.output["properties"]
        jobs.append(Np1)

        calculator.charge = baseCharge - 1
        Nm1 = calculator.make(structure=structure)
        Nm1_energy = Nm1.output["properties"]
        jobs.append(Nm1)

        calculate = CalculateGlobalConceptualDFTProperties().make(
            N_energy, Np1_energy, Nm1_energy
        )
        jobs.append(calculate)

        flow = Flow(jobs, name="CDFT Calculation")
        return Response(replace=flow)
