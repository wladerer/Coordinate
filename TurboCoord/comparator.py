from xyz2graph import MolGraph, to_networkx_graph, to_plotly_figure
import networkx as nx
import matplotlib.pyplot as plt
import networkx.algorithms.isomorphism as iso
import itertools
from itertools import permutations
import os


def graph(file):
    mg = MolGraph()
    mg.read_xyz(file)
    G = to_networkx_graph(mg)
    return G

def compareGraphs(G1,G2):

    return nx.is_isomorphic(G1, G2)


directory = '/home/wladerer/github/Coordinate/utils/smaller_sample'

graphs = []
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    graphs.append(graph(f))


unique_combinations = []
 
# Getting all permutations of list_1
# with length of list_2
permut = itertools.permutations(graphs, len(graphs))
 
# zip() is called to pair each permutation
# and shorter list element into combination
for comb in permut:
    zipped = zip(comb, graphs)
    unique_combinations.append(list(zipped))


is_iso = 0
for tup in unique_combinations:
    tup = tup[0]
    if compareGraphs(tup[0],tup[1]):
        is_iso += 1

