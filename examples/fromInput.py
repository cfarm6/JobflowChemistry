import json
import networkx as nx
import matplotlib.pyplot as plt
import jobflowchemistry
from jobflow import Flow
from jobflow.managers.fireworks import flow_to_workflow
from fireworks import LaunchPad
from icecream import ic
import uuid

with open("input.json", "r") as f:
    json_data = json.load(f)

del json_data["viewport"]

input_nodes = list(filter(lambda x: x["type"] != "jobflowNode", json_data["nodes"]))
input_ids = list(map(lambda x: x["id"], input_nodes))

jobflow_nodes = list(filter(lambda x: x["type"] == "jobflowNode", json_data["nodes"]))
jobflow_ids = list(map(lambda x: x["id"], jobflow_nodes))
G = nx.MultiDiGraph()


def initMaker(node):
    # Get maker with initial values
    maker = getattr(
        __import__(f"jobflowchemistry.{node["data"]["parent"]}"), node["data"]["type"]
    )
    # Set values for maker
    val = map(
        lambda x: {x: node["data"]["settings"][x]["value"]}, node["data"]["settings"]
    )
    settings = {}
    list(map(lambda x: settings.update(x), val))
    maker.from_dict(settings)
    ic(maker)
    # Check node settings
    return maker()


# Create input nodes
for node in input_nodes:
    G.add_node(node["id"], data=node["data"])

# Create jobflow nodes and populate settings
for node in jobflow_nodes:
    nodeMaker = initMaker(node)
    G.add_node(node["id"], node=nodeMaker)

# Add the edges to the graph
for edge in json_data["edges"]:
    G.add_edge(
        edge["source"],
        edge["target"],
        sourceHandle=edge["sourceHandle"],
        targetHandle=edge["targetHandle"],
    )

def makeNode(nodeId):
    # Skip if the node is in input_ids
    if nodeId in input_ids: return
    # Skip if the node was already made
    if "job" in G.nodes[nodeId]: return
    input_dict = {}
    # Get the predecessors
    for pred in G.predecessors(nodeId):
        edges = G.get_edge_data(pred, nodeId)
        for e in edges:
            # Check if it is an input node and add the inputs
            if pred in input_ids:
                input_dict[edges[e]["targetHandle"]] = G.nodes[pred]["data"][
                    edges[e]["sourceHandle"]
                ]
                continue

            # If the Node is being used connect the node output instead of making the node
            elif edges[e]["sourceHandle"] == "Node":
                input_dict[edges[e]['targetHandle']] = G.nodes[pred]['node']
                continue

            # Check if the job has been made yet and make it if it has not been made
            elif "job" not in G.nodes[pred]:
                makeNode(pred)

            input_dict[edges[e]["targetHandle"]] = G.nodes[pred]["job"].output[
                edges[e]["sourceHandle"]
            ]

    # Conditions for not making the job:
    ## 1. No input connections AND
    ## 2. The only output connection is connected to 'Node'
    if len(G.in_edges(nodeId)) == 0: 
        if all(map(lambda x: all(map(lambda y: G[nodeId][x][y]['sourceHandle']=='Node', G[nodeId][x])), G[nodeId])):
            return 
        
    G.nodes[nodeId]["job"] = G.nodes[nodeId]["node"].make(**input_dict)
    ic(G.nodes[nodeId]["job"])
    return

# 1. Make jobs from the nodes
for node in jobflow_nodes:
    makeNode(node["id"])

jobs = [G.nodes[x]["job"] for x in jobflow_ids if "job" in G.nodes[x]]
ic(len(jobs))
# 2. Get the jobs for the flow
flow = Flow(jobs, name = str(uuid.uuid4()))
fw = flow_to_workflow(flow)
lpad = LaunchPad.auto_load()
lpad.add_wf(fw)
