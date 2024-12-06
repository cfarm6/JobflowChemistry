from ..PropertyCalculator import (
    CalculateGlobalConceptualDFTProperties,
    CollisionCrossSectionCalculator,
)
from ..EnergyCalculation import EnergyCalculation
from ..outputs import Settings, Properties
from ..Structure import Structure
from ..ConformerGeneration import ConformerGeneration
from ..GeometryOptimization import GeometryOptimization
from ..Utilities import BoltzmannWeighting
from rdkit.Chem import rdmolops, rdchem
from dataclasses import dataclass, field
from jobflow import Response, job, Flow, Maker
from typing import Literal

@dataclass
class PropertyWorkflows(Maker):
    name: str = "Property Workflows"
    @job 
    def make():
        raise NotImplementedError


@dataclass
class ConceptualDFTWorkflow(PropertyWorkflows):
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


@dataclass
class CollisionCrossSectionWorkflow(PropertyWorkflows):
    name: str = "Collision Cross Section Workflow"
    temperature: float = 300.0
    property_type: Literal["global"] = "global"
    property: Literal["CCS [A^2]"] = "CCS [A^2]"
    energy_name: str = "Total Energy [eV]"

    @job(files="files", settings="settings", properties="properties")
    def make(
        self,
        structure: Structure,
        calculator_conformer: ConformerGeneration,
        calculator_optimization: GeometryOptimization,
        calculator_energy: EnergyCalculation,
        calculator_ccs: CollisionCrossSectionCalculator,
    ):
        if type(structure) is list:
            jobs = [
                self.make(
                    s, calculator_conformer, calculator_optimization,calculator_energy, calculator_ccs
                )
                for s in structure
            ]
            return Response(
                output={
                    "structure": [x.output["structure"] for x in jobs],
                    "settings": Settings({}),
                    "properties": [x.output["properties"] for x in jobs],
                },
                addition=jobs,
            )
        if structure.GetNumConformers() > 1:
            jobs = []
            for confId in range(structure.GetNumConformers()):
                s = Structure(rdchem.Mol(structure, confId=confId))
                jobs.append(
                    self.make(
                        s, calculator_conformer, calculator_optimization,calculator_energy, calculator_ccs
                    )
                )
            return Response(
                output={
                    "structure": [x.output["structure"] for x in jobs],
                    "settings": Settings({}),
                    "properties": [x.output["properties"] for x in jobs],
                },
                addition=jobs,
            )
        boltzmann_calculator = BoltzmannWeighting(
            temperature=self.temperature,
            property=self.property,
            property_type=self.property_type,
            energy_name=self.energy_name,
        )

        conformer_job = calculator_conformer.make(structure=structure)
        optimization_job = calculator_optimization.make(
            structure=conformer_job.output["structure"]
        )
        energy_job = calculator_energy.make(structure=optimization_job.output["structure"])
        ccs_job = calculator_ccs.make(structure=energy_job.output["structure"])
        boltzmann_job = boltzmann_calculator.make(structure=ccs_job.output["structure"])
        return Response(
            addition=Flow([conformer_job, optimization_job, energy_job, ccs_job, boltzmann_job])
        )
