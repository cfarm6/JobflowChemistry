from jobflowchemistry.StructureInput import PubChemInput
from jobflowchemistry.StructureGeneration import (
    RDKitGeneration,
)
from jobflowchemistry.GeometryOptimization import AimNet2Optimization
from jobflowchemistry.EnergyCalculation import xTBEnergyCalculation, PTBCalculation

from jobflowchemistry.HessianCalculation import xTBHessian
from jobflowchemistry.SpectraCalculator import PTBRamanSpectra
#
from jobflow.managers.fireworks import flow_to_workflow
from fireworks import LaunchPad
from jobflow import Flow, job
import pandas as pd
from icecream import ic

# Node Settings
input_maker = PubChemInput()
structure_maker = RDKitGeneration()
aimnet_maker = AimNet2Optimization()
xtb_maker = xTBHessian(solvent="dmf", cosmo=True)
raman_maker = PTBRamanSpectra()
# Make the workflow
input_job = input_maker.make(input=10784527)
gen_job = structure_maker.make(structure=input_job.output["structure"])
opt_job = aimnet_maker.make(structure=gen_job.output["structure"])
raman_job = raman_maker.make(structure=opt_job.output["structure"])
flow = Flow(
    [
        input_job,
        gen_job,
        opt_job,
        raman_job,
    ],
    name=str(10784527),
)
fw = flow_to_workflow(flow)
lpad = LaunchPad.auto_load()
lpad.add_wf(fw)
