"""This file contains function which help in parsing GPML file."""
import networkx as nx
from bs4 import BeautifulSoup
from itertools import permutations, product


#oggpnosn
#hkhr

def get_node_information(data_node):
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



class PathwayGraph:

    def __init__(self, pathway_filepath):
        """
        It creates a pathway graph by parsing GPML file.
        #Param:
            pathway_filepath: filepath of the GPML file.
        """
        self.graph = self._get_graph(pathway_filepath)
        # Pre-compute leaf and root nodes for faster matching.
        self.leaf_nodes = [node_id for node_id in self.graph.nodes \
                           if self.graph.in_degree(node_id) != 0 and \
                           self.graph.out_degree(node_id) == 0]
        self.root_nodes = [node_id for node_id in self.graph.nodes \
                           if self.graph.in_degree(node_id) == 0 and \
                           self.graph.out_degree(node_id) != 0]



    def _find_target_in_graph(self, ensemble_ids):
        nodes = set()
        # Only check root or leaf nodes since everything in between is
        # just a transport mechanism.
        for node in self.leaf_nodes + self.root_nodes:
            if self.graph.nodes[node]["attr_dict"]["type"] == "Complex":
                # In case of protein complex. Check if any of the esemble id
                # lies in the protein complex.
                database_ids = self.graph.nodes[node]["attr_dict"]["member_nodes"]
                for ensemble_id in ensemble_ids:
                    if ensemble_id in database_ids:
                        nodes.add(node)
            else:
                # In case of a non complex node, just check if database_id is
                # any of the ensemble_ids
                database_id = self.graph.nodes[node]["attr_dict"].get("database_id", "<unk>")
                if database_id in ensemble_ids:
                    nodes.add(node)

        return list(nodes)

    def get_shortest_path(self, gene_1, gene_2):
        """
        It finds the shortest path between the two genes if there exist any. In case
        there exist no path between two genes it returns empty array. For the cases
        where the path exist, it returns the path as well as path length.
        """
        nodes_1 = self._find_target_in_graph(gene_1)
        nodes_2 = self._find_target_in_graph(gene_2)

        path_information = []
        for node_1, node_2 in product(nodes_1, nodes_2):
            if node_1 ==  node_2:
                # It makes no sense to find the shortest path for this case.
                continue
            try:
                path = nx.shortest_path(self.graph, source=node_1,
                                        target=node_2)
                # Remove all anchors from the path.
                path = [self.graph.node[node_id]["attr_dict"]["text"] for node_id in path \
                        if self.graph.node[node_id].get("attr_dict", {}).get("type", "") not in ["anchor", ""]]
                path_information.append({
                    "path": path,
                    "path_length": len(path)-1
                })
            except nx.exception.NetworkXNoPath:
                continue

        return path_information

    def get_shortest_path_for_list(self, database_ids, pathway_name):
        """
        It returns shortest path between all the pairs of gene's in the list.
        :param database_ids: a map of gene name to it's corresponding database ids like ensemble id or chembl id.
        :param pathway_name: the name of the pathway that goes into constraints.
        :return: path information for all pair of gene's in the gene_list
        """
        # Find the nodes of interest.
        nodes_of_interest = set()
        graph_map = {} # map of graph's id to gene_names.
        for gene_name in database_ids:
            gene_nodes = self._find_target_in_graph(database_ids[gene_name])
            nodes_of_interest = nodes_of_interest.union(set(gene_nodes))
            for gene_node_id in gene_nodes:
                node_mapping = graph_map.get(gene_node_id, [])
                node_mapping.append(gene_name)
                graph_map[gene_node_id] = node_mapping

        path_constraints = []
        for gene_node in nodes_of_interest:
            # Using nodes of interest as source, compute single source shortest path algorithm.
            path_map = nx.single_source_shortest_path(self.graph, source=gene_node)
            # for each node of interest retain the information in the same format as before.
            for target_node in nodes_of_interest:
                if len(path_map.get(target_node, [])) > 1:
                    path_to_target = path_map[target_node]
                    # Remove all anchors from the path.
                    path = [self.graph.node[node_id]["attr_dict"]["text"] for node_id in path_to_target \
                            if self.graph.node[node_id].get("attr_dict", {}).get("type", "") not in ["anchor", ""]]

                    for gene_source_name, gene_sink_name in product(graph_map[gene_node], graph_map[target_node]):
                        path_constraints.append((gene_source_name, gene_sink_name, {
                            "type": "pathway_path",
                            "path_information": [{"path": path,
                                                 "path_length": len(path) - 1}],
                            "pathway_name": pathway_name,
                            "meta_info": self.pathway_name,
                            "color_id":pathway_name,
                        }))
        return path_constraints


    def get_immediate_nodes(self, gene_1,hops=1):
        """
        It finds immediate nodes of gene_1
        """
        nodes_1 = self._find_target_in_graph(gene_1)
        path_information = []
        for node_1 in nodes_1:
            try:
                neighbor_paths =nx.single_source_shortest_path(self.graph, node_1, cutoff=hops*2)
                neighbor_paths = [ j for i,j in neighbor_paths.items()]
                processed_neighbour_paths = []
                for path in neighbor_paths:
                    path = [self.graph.node[node_id]["attr_dict"]["text"] for node_id in path \
                            if self.graph.node[node_id].get("attr_dict", {}).get("type", "") not in ["anchor", ""]]
                    if len(path) <= hops:
                        # for neighbor_node in path[-1].split('|'):
                        #     processed_neighbour_paths.append(path + [neighbor_node])
                        processed_neighbour_paths.append(path )
                processed_neighbour_paths.append([node_1])
                neighbor_paths = processed_neighbour_paths
                path_information.extend([{
                    "path": path,
                    "path_length": len(path)-1
                } for path in neighbor_paths])
            except nx.exception.NetworkXNoPath:
                continue
        return path_information


if __name__ == "__main__":
    graph = get_networkx_graph("../Data/Hs_p75_NTR_receptor-mediated_signalling_WP4443_101574.gpml")

    print(graph.number_of_nodes())
    print(graph.number_of_edges())

    print(graph.node["ba152"])

    print(graph.node["bb7b1"])

    nx.write_gexf(graph, "sample_graph.gexf")
    from collections import Counter
    print(Counter([graph.node[node_id]["attr_dict"]["type"] for node_id in graph.nodes]))


