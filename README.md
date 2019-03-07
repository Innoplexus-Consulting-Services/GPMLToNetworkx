#oggpnosn 
#hkhr 




GPML Pathway to networkx convertor
========

GPMLToNetworkx is a Python package for parsing GPML pathway files into a 
networkx graph for downstream processing.


Install
-------

Install the latest version of GPMLToNetworkx::

    $ pip install .

Simple example
--------------

Read the GPML file into a networkx graph:


    >>> from GPMLParser import get_networkx_graph
    >>> graph = get_networkx_graph("./Data/Hs_p75_NTR_receptor-mediated_signalling_WP4443_101574.gpml")
    >>> graph.number_of_nodes()
    244
    >>> graph.number_of_edges()
    341
    >>> print(graph.node["ba152"])
    {'attr_dict': {'database_name': 'Uniprot-TrEMBL', 'database_id': 'Q13501', 'text': 'SQSTM1', 'type': 'Protein', 'groupref': None}}
    >>> print(graph.node["bb7b1"])
    {'attr_dict': {'member_nodes': ['P63000', 'CHEBI:15996'], 'type': 'Complex', 'text': 'GGC-PalmC-RAC1 |GTP '}}
    
Bugs
----

License
-------

