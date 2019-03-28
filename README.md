# GPML Pathway to networkx convertor

GPMLToNetworkx is a Python package for parsing GPML pathway files into a 
networkx graph for downstream processing.


## Install


    $ pip install -r requirements.txt

## Simple example


Read the GPML file into a networkx graph:
```python
from GPMLParser import get_networkx_graph
graph = get_networkx_graph("./Data/Hs_p75_NTR_receptor-mediated_signalling_WP4443_101574.gpml")
```    
Once you have networkx graph object, you can use all standard networkx operations on it.

    >>> graph.number_of_nodes()
    244
    >>> graph.number_of_edges()
    341
    >>> print(graph.node["ba152"])
    {'attr_dict': {'database_name': 'Uniprot-TrEMBL', 'database_id': 'Q13501', 'text': 'SQSTM1', 'type': 'Protein', 'groupref': None}}
    >>> print(graph.node["bb7b1"])
    {'attr_dict': {'member_nodes': ['P63000', 'CHEBI:15996'], 'type': 'Complex', 'text': 'GGC-PalmC-RAC1 |GTP '}}
 
## Node Semantics

Nodes can be of three types:
1. **Protein/Gene/Metabolite/Pathway** - These nodes refer to only one entity
2. **Complex** - These nodes refer to a collection of entities that have formed complex.
3. **Anchor** - These nodes provide connectivity. They must not be taken into consideration in Algorithm. Since connections
in  pathway aren't edges in graph theoretic sense, we converted the connnection into 
group anchor nodes that are connected to one another. You can see the shortest path 
computation to see this point in action.

## Edge Semantics

There are two types of connection:
1. **Non-Synthetic**: This connection was originally present in the pathway.
2. **Synthetic**: This connection was introduced to make sure anchor nodes in a connection are connected to one another.

## Shortest path example

If you want to find shortest path length between two non-anchor nodes, you must ignore all 
the anchor nodes. Here's how you can do it:
```python
import networkx as nx 

path = nx.shortest_path(graph, source=source_node, target=target_node)
refined_path = [graph.node[node_id]["attr_dict"]["text"] for node_id in path if graph.node[node_id]["attr_dict"]["type"] != "anchor"]
path_length = len(refined_path)-1
```    

## Bugs
You can raise issues in github issues.

## Authors

Tanay Gahlot - tanay.gahlot@innoplexus.com


