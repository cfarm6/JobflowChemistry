from jobflowchemistry.StructureInput import PubChemInput
from jobflowchemistry.StructureGeneration import (
    RDKitGeneration,
    CRESTProtonation,
    CRESTDeprotonation,
    aISSDocking,
)
from jobflowchemistry.EnergyCalculation import (
    AimNet2EnergyCalculation,
    ORCAEnergyCalculation,
)
from jobflowchemistry.GeometryOptimization import (
    AimNet2Optimization,
    ORCAOptimization,
)
from jobflowchemistry.PropertyCalculator import (
    QupKakePrediction,
    ConceptualDFTWorkflow,
)
from jobflowchemistry.ConformerGeneration import (
    rdKitConformers,
    Auto3DConformers,
    CRESTConformers,
)

#
from jobflow.managers.fireworks import flow_to_workflow
from fireworks import LaunchPad
from jobflow import Flow, job
import pandas as pd
from icecream import ic

# Node Settings
input_maker = PubChemInput()
structure_maker = RDKitGeneration()
orcaenergy_maker = ORCAEnergyCalculation(
    parallel="PAL16", 
    scf_conv="SLOPPYSCF",
    executable="/home/carson/orca_6_0_1/orca"
)
aimnet2energy_maker = AimNet2EnergyCalculation()
aimnet2opt_maker = AimNet2Optimization()
orcaopt_maker = ORCAOptimization(
    parallel="PAL16",
    geom_conv="LOOSEOPT",
    scf_conv="SLOPPYSCF",
    executable="/home/carson/orca_6_0_1/orca",
)
cdft_calculation = ConceptualDFTWorkflow()
qupkake_maker = QupKakePrediction(mp=False)
protonation_maker = CRESTProtonation(ion='Na+', ion_charge=1, ewin = 10.0, threads=20)
deprotonation_maker = CRESTDeprotonation(ewin = 30.0)
aISSDocking_maker = aISSDocking(fast = True, ensemble=False)
rdKitConf_maker = rdKitConformers(numConfs=2)
auto3d_maker = Auto3DConformers(window = 10.0)
crestconf_maker = CRESTConformers(energy_method='gfn2', dynamics_method='gfnff', threads=20, ewin = 3.0)
# Build workflows
# df = pd.read_csv("examples/cdft_fluorophenols.csv")
# for index, row in df.iterrows():
# if index != 0: continue
input_job = input_maker.make(input=75921)
gen_job = structure_maker.make(structure=input_job.output["structure"])
aimnet2opt_job = aimnet2opt_maker.make(structure=gen_job.output["structure"])
rdKitConf_job = rdKitConf_maker.make(structure=aimnet2opt_job.output["structure"])
auto3d_job = auto3d_maker.make(structure=aimnet2opt_job.output["structure"])
crestconf_job = crestconf_maker.make(structure=aimnet2opt_job.output["structure"])
# orcaenergy_job = orcaenergy_maker.make(structure=aimnet2opt_job.output['structure'])
# orcaopt_job = orcaopt_maker.make(structure=aimnet2opt_job.output["structure"])
# cdft_job = cdft_calculation.make(
#     structure=aimnet2opt_job.output["structure"],
#     neutral_properties=aimnet2opt_job.output["properties"],
#     calculator=aimnet2energy_maker,
# )
# deprotonate_job = deprotonation_maker.make(
#     structure=aimnet2opt_job.output["structure"],
# )
# protonate_job = protonation_maker.make(
#     structure=deprotonate_job.output["structure"],
# )
# docking_job = aISSDocking_maker.make(
#     structure_1=deprotonate_job.output["structure"],
#     structure_2=protonate_job.output["structure"],
# )
# qupkake_job = qupkake_maker.make(structure=aimnet2opt_job.output["structure"])
#
flow = Flow(
    [
        input_job,
        gen_job,
        # orcaenergy_job,
        aimnet2opt_job,
        crestconf_job,
        # orcaopt_job,
        # cdft_job,
        # qupkake_job,
        # protonate_job,
        # deprotonate_job,
        # docking_job,
        rdKitConf_job,
        # aimnet2energy_job,
        # aimnet2opt1_job,
        auto3d_job,
        # aimnet2energy2_job,
    ],
    name=str(2776882),
)
fw = flow_to_workflow(flow)
lpad = LaunchPad.auto_load()
lpad.add_wf(fw)
