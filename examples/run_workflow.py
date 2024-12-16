from jobflowchemistry.StructureInput import PubChemInput
from jobflowchemistry.StructureGeneration import (
    RDKitGeneration,
    CRESTDeprotonation,
)
from jobflowchemistry.EnergyCalculation import (
    AimNet2EnergyCalculation,
    TBLiteEnergyCalculation
)
from jobflowchemistry.GeometryOptimization import (
    AimNet2Optimization,
)
from jobflowchemistry.PropertyCalculator import (
    MobcalCollisionCrossSection,
)
from jobflowchemistry.ConformerGeneration import (
    Auto3DConformers,
)
from jobflowchemistry.PropertyWorkflows import CollisionCrossSectionWorkflow

#
from jobflow.managers.fireworks import flow_to_workflow
from fireworks import LaunchPad
from jobflow import Flow, job
import pandas as pd
from icecream import ic

# Node Settings
input_maker = PubChemInput()
structure_maker = RDKitGeneration()
aimnet2opt_maker = AimNet2Optimization()
aimnet2energy_maker = TBLiteEnergyCalculation()
deprotonation_maker = CRESTDeprotonation(ewin=30.0, checkStructure=False)
auto3d_maker = Auto3DConformers(k=2)
mobcal_maker = MobcalCollisionCrossSection(
    charge_type="Mulliken Partial Charges [e]", num_threads=5, temperature=300.0
)
ccs_workflow = CollisionCrossSectionWorkflow(temperature=300.0)
# Make the workflow
input_job = input_maker.make(input=10784527)
gen_job = structure_maker.make(structure=input_job.output["structure"])
tblite_job = aimnet2energy_maker.make(structure=input_job.output["structure"])
aimnet2opt_job = aimnet2opt_maker.make(structure=gen_job.output["structure"])
deprotonate_job = deprotonation_maker.make(
    structure=aimnet2opt_job.output["structure"],
)
ccs_job = ccs_workflow.make(
    structure=deprotonate_job.output["structure"],
    conformer_calculator=auto3d_maker,
    optimization_calculator=aimnet2opt_maker,
    energy_calculator=aimnet2energy_maker,
    ccs_calculator=mobcal_maker,
)

flow = Flow(
    [
        input_job,
        gen_job,
        tblite_job,
        aimnet2opt_job,
        # deprotonate_job,
        # ccs_job,
    ],
    name=str(75921),
)
fw = flow_to_workflow(flow)
lpad = LaunchPad.auto_load()
lpad.add_wf(fw)
