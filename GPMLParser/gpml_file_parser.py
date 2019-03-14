"""This file contains function which help in parsing GPML file."""
import networkx as nx
from bs4 import BeautifulSoup
from itertools import permutations, product


#oggpnosn
#hkhr

def get_node_information(data_node):
    """
    It parses nodes in GPML xml file.
    :param data_node: XML node of the type "data"
    :return: node id and structured information of the form:
    {
        "database_name": database_name,
        "database_id": database_id,
        "text": text,
        "type": type,
        "groupref": groupref
    }
    """
    # get database and its corresponding unique identifier.
    xref = data_node.find("xref")
    database_name = xref.attrs["database"]
    database_id = xref.attrs["id"]

    # get the graph id.
    id = data_node.attrs.get("graphid", None)
    text = data_node.attrs["textlabel"]
    type = data_node.attrs.get("type", "unknown")
    groupref = data_node.attrs.get("groupref", None)

    # form a document.
    return (id, {
        "database_name": database_name,
        "database_id": database_id,
        "text": text,
        "type": type,
        "groupref": groupref
    })


def get_edge_information(interaction):
    """
    It parses the interaction node in GPML file.
    :param interaction: the XML node whose type is <interaction>
    :return: structured information about the interaction. It's of the form:
    {
        "source_ids": source_ids,
        "sink_id": sink_id,
        "edge_type": edge_type,
        "anchors": anchors
    }
    """
    graphics_element = interaction.find("graphics")
    points = graphics_element.find_all("point")

    source_ids = []
    for point in points[:-1]:
        source_id = point.attrs.get("graphref", False)
        if source_id:
            source_ids.append(source_id)

    sink_id = points[-1].attrs.get("graphref", None)
    edge_type = points[-1].attrs.get("arrowhead", "unknown")

    # finding all the anchor ids.
    anchors = graphics_element.find_all("anchor")
    anchors = [anchor.attrs["graphid"] for anchor in anchors if "graphid" in anchor.attrs]

    return {
        "source_ids": source_ids,
        "sink_id": sink_id,
        "edge_type": edge_type,
        "anchors": anchors
    }


def get_group_information(group_node, data_nodes):
    """
    It parses the complex in the GPML file.
    :param group_node: the nodes in the group
    :param data_nodes: all the data nodes present in the GPML file.
    :return: group id and structured information about the group.
    {"member_nodes": member_nodes, "type": "Complex",
     "text": text}
    """
    graph_id = group_node.attrs.get("graphid", None)
    group_ref = group_node.attrs["groupid"]
    attribute = group_node.find("attribute")

    # Find all the members of the complex.
    member_nodes = []
    text_labels = []
    for data_node in data_nodes:
        if data_node.attrs.get("groupref", None) == group_ref:
            text = data_node.attrs["textlabel"]
            xref = data_node.find("xref")
            ensemble_id = xref.attrs["id"]
            if ensemble_id:
                member_nodes.append(ensemble_id)
                text_labels.append(text)

    # In this wicked scenario, the graph_id doesnt refers to the real graph id.
    # we would have to find out the real graph id by matching
    if attribute:
        complex_id = attribute.attrs["value"]
        for data_node in data_nodes:
            xref = data_node.find("xref")
            db_id = xref.attrs["id"]
            if db_id == complex_id:
                graph_id = data_node.attrs.get("graphid", None)

    return graph_id, {"member_nodes": member_nodes, "type": "Complex",
                       "text": "|".join(text_labels)}


def get_networkx_graph(pathway_filepath):
    """it converts the GPML format to networkx graph"""
    graph = nx.DiGraph()

    wiki_file = BeautifulSoup(open(pathway_filepath), "html.parser")

    data_nodes = wiki_file.find_all("datanode")
    for data_node in data_nodes:
        node_id, node_information = get_node_information(data_node)
        # node_id can be empty at times.
        if node_id and not node_information["groupref"]:
            graph.add_node(node_id, attr_dict=node_information)

    group_nodes = wiki_file.find_all("group")
    for group_node in group_nodes:
        group_id, group_info = get_group_information(group_node, data_nodes)
        if group_id:
            graph.add_node(group_id, attr_dict=group_info)

    interactions = wiki_file.find_all("interaction")
    for interaction in interactions:
        edge_information = get_edge_information(interaction)

        # Creating an edge between every source and sink node.
        sink_id = edge_information["sink_id"]
        edge_type = edge_information["edge_type"]
        for source_id in edge_information["source_ids"]:
            # In case the source or sink is anchor node, it wouldn't exist
            # in the graph. Therefore before we add an edge, we have to
            # introduce anchor.
            if source_id not in graph:
                graph.add_node(source_id, attr_dict={"type": "anchor"})
            if sink_id not in graph:
                graph.add_node(sink_id, attr_dict={"type": "anchor"})
            graph.add_edge(source_id, sink_id, attr_dict={"edge_type": edge_type})

        # Ensure all anchors are in the graph.
        for anchor_id in edge_information["anchors"]:
            if anchor_id not in graph:
                graph.add_node(anchor_id, attr_dict={"type": "anchor"})
        # Connecting all the anchors with each other.
        for anchor_1, anchor_2 in permutations(edge_information["anchors"], 2):
            graph.add_edge(anchor_1, anchor_2)
    return graph
