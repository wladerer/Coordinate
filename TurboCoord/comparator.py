import networkx as nx
import matplotlib.pyplot as plt
import networkx.algorithms.isomorphism as iso
import itertools
from itertools import permutations
import os
from TurboCoord import CoordinationComplex, generateStructures

class MolecularGraph:
    """
    Represent the xyz coordinates of a molecule in terms of a graph using networkx
    """

    def __init__(self, CoordinationComplex):
        super().__init__(CoordinationComplex)
        self.graph = nx.Graph()
        self.graph.add_nodes_from(self.complex_atoms)
        self.graph.add_edges_from(self.get_edges())
        self.graph = self.graph.to_undirected()
        self.graph = nx.convert_node_labels_to_integers(self.graph, first_label=0, ordering='default', label_attribute=None)

    def get_edges(self):
        """
        Get the edges of the graph
        """
        edges = []
        for i in range(len(self.complex_atoms)):
            for j in range(len(self.complex_atoms)):
                if self.distance_matrix[i][j] < 2.0:
                    edges.append((self.complex_atoms[i], self.complex_atoms[j]))
        return edges

    def get_graph(self):
        return self.graph

    def get_isomorphisms(self, other):
        """
        Get the isomorphisms between two graphs
        """
        GM = iso.GraphMatcher(self.graph, other.graph)
        return GM.is_isomorphic()
    
