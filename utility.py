import numpy as np
from igraph import *
import random
import copy


def get_graph(edge_list):
    
    
    g = Graph.Read_Ncol(edge_list, names=True, weights="if_present", directed=False)
    
    return g

def get_hash(graph):
    neighbor_dict = {}
    indices = graph.vs.indices
    for index in indices:
        neighbor = graph.neighborhood(index)
        neighbor.remove(index)
        neighbor_dict[index] = neighbor
    return indices, neighbor_dict
