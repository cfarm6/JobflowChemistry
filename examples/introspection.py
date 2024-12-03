import importlib
import inspect
import networkx as nx
import matplotlib.pyplot as plt
from jobflow import Maker
import json
from dataclasses import fields, MISSING
from icecream import ic
from typing import get_type_hints, Literal, get_args

def get_inputs(node):
    inputs = []
    type_hints = get_type_hints(node.make)
    for type_hint in type_hints:
        inputs.append({"name": type_hint, "type": type_hints[type_hint].__name__})
    return inputs


def get_settings(node):
    settings = fields(node)
    settings_dict = {}
    for setting in settings:
        # if setting.name == "name":
        #     continue
        _setting = {
            "type": str(setting.type.__name__),
        }
        if setting.type.__name__ == "Literal":
            _setting["options"] = get_args(setting.type)
        else:
            _setting["options"] = None
        if setting.default is not MISSING:
            _setting["value"] = setting.default
        else:
            _setting["value"] = None
        settings_dict[setting.name] = _setting
    if settings_dict == {}:
        return None
    name = settings_dict['name']['value']
    del settings_dict['name']
    options = {
        "settings": settings_dict,
        "inputs": get_inputs(node),
        "name": name,
        "type": node.__name__
    }
    return options


def build_class_tree(module_name):
    # Dynamically import the module
    module = importlib.import_module(module_name)

    # Initialize a directed graph
    G = nx.DiGraph()

    # Package Classes
    classes = [
        (name, cls)
        for name, cls in inspect.getmembers(module, inspect.isclass)
        if issubclass(cls, Maker) and cls is not Maker
    ]
    # Traverse through all items in the module
    for name, obj in classes:
        # Add the class as a node in the graph
        get_inputs(obj)
        G.add_node(obj.__name__, data = get_settings(obj))
        # Loop through the class's base classes
        for base in obj.__bases__:
            if base in classes: continue
            # Add an edge from the base class to the current class
            G.add_edge(base.__name__, obj.__name__)

    return G


def get_highest_level_nodes(G):
    # Find all nodes with in-degree of 0 (no parents)
    highest_level_nodes = [node for node in G.nodes if G.in_degree(node) == 0]
    return highest_level_nodes


def remove_non_maker_highest_level_nodes(G):
    # Find all nodes with in-degree of 0 (no parents)
    highest_level_nodes = [node for node in G.nodes if G.in_degree(node) == 0]

    # Filter out nodes that do not contain 'Maker' in their name
    for node in highest_level_nodes:
        if "Maker" not in node:
            G.remove_node(node)

    return G


def plot_class_tree(G):
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G)
    nx.draw(
        G,
        pos,
        with_labels=True,
        node_size=2000,
        node_color="lightblue",
        font_size=10,
        font_weight="bold",
        edge_color="gray",
    )
    plt.title("Class Inheritance Tree")
    plt.show()


# Example Usage
module_name = "jobflowchemistry"  # replace with your package's name
class_tree = build_class_tree(module_name)


# Example usage
highest_level_nodes = get_highest_level_nodes(class_tree)
print("Highest level nodes (root classes):", highest_level_nodes)


# Example usage
class_tree_filtered = remove_non_maker_highest_level_nodes(class_tree)
print("Remaining nodes after filtering:", class_tree_filtered.nodes())
highest_level_nodes = get_highest_level_nodes(class_tree_filtered)
print("Highest level nodes (root classes):", highest_level_nodes)
plot_class_tree(class_tree_filtered)
with open("output.json", "w") as f:
    f.write(
        json.dumps(nx.node_link_data(class_tree_filtered, name="name", key="settings"))
    )
