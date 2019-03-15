"""it ensures that the GPML parser works well."""
#oggpnosn
#hkhr
from GPMLParser import get_networkx_graph
import networkx as nx
import pytest


def get_shortest_path_length(graph, source_node, target_node):
    path = nx.shortest_path(graph, source=source_node, target=target_node)
    refined_path = [graph.node[node_id]["attr_dict"]["text"] for node_id in path if
                    graph.node[node_id]["attr_dict"]["type"] != "anchor"]
    path_length = len(refined_path) - 1
    return path_length


def get_node(graph, identifier):
    for node in graph.nodes:
        if graph.node[node]["attr_dict"].get("database_id", None) == identifier:
            return node
    return None


class TestGPMLParser:


    def test_shortest_path(self):
        source_node = get_node(graph, "Q9Y4K3")
        target_node = get_node(graph, "P08138")
        with pytest.raises(nx.exception.NetworkXNoPath):
            get_shortest_path_length(graph, source_node, target_node)

        source_node = get_node(graph, "P45983")
        target_node = get_node(graph, "Q92934")
        with pytest.raises(nx.exception.NetworkXNoPath):
            get_shortest_path_length(graph, source_node, target_node)


    def test_node_existence(self):
        assert get_node(graph, "Q9Y4K3") is not None
        assert get_node(graph, "P41743") is not None
        assert get_node(graph, "Q13352") is not None
        assert get_node(graph, "O43521") is not None
        assert get_node(graph, "P45983") is not None



graph = get_networkx_graph("../Data/Hs_p75_NTR_receptor-mediated_signalling_WP4443_101574.gpml")