from utility import get_graph, get_hash
import networkx as nx
import pickle
import sys


def lcc_extraction(edge_list):
    "extract the largest connected component"
    with open(edge_list, 'r') as fi:
        lines = fi.readlines()
    del lines[0]
    H = nx.Graph()
    edges = []
    for line in lines:
        splitted_line = line.split()
        edge = (int(splitted_line[0]), int(splitted_line[1]))
        edges.append(edge)
    
    H.add_edges_from(edges)
    G = H.to_undirected()
    largest = max(nx.connected_component_subgraphs(G), key = len)
    nx.write_edgelist(largest, edge_list + '_lcc')
    
def trim_edges(edge_list_lcc):
    "trim edges"
    with open(edge_list_lcc, 'r') as fi:
        lines = fi.readlines()
    node_dict = {}
    i = 1
    for line in lines:
        splitted_line = line.split()
        for word in splitted_line:
            if word not in node_dict:
                if word != '{}':
                    node_dict[word] = i
                    i += 1
    new_lines = []
    for line in lines:
        new_line = ''
        splitted_line = line.split()
        new_line += str(node_dict[splitted_line[0]]) + '\t' + str(node_dict[splitted_line[1]]) + '\n'
        new_lines.append(new_line)
    with open(edge_list_lcc + '_trimmed', 'w') as fo:
        fo.writelines(new_lines)
        
def delete_self_loop(edgelist, outputfilename):
    "delete self loops"
    g = get_graph(edgelist)
    length = len(g.es)
    
    while True :
        selfloop = False;
        for e in g.es:
            if g.is_loop(e):
                selfloop = True;
                g.delete_edges(e)
        if selfloop == False : 
            break;    
    
    g.write_edgelist(outputfilename)

if __name__ == "__main__":
    edgelist = sys.argv[1]
    print('extracting the largest connected component')
    lcc_extraction(edgelist)
    print('trimming edges')
    trim_edges(edgelist + '_lcc')
    print('deleting self loops')
    delete_self_loop(edgelist + "_lcc", edgelist + "_final")
    print('getting the graph to make hash file of edgelist')
    g = get_graph(edgelist + "_final")
    
    print('hashing the edgelist')
    graph_hash = get_hash(g)
    print('dumping the edgelist')
    hashfile = edgelist + '_hash.p'
    pickle.dump(graph_hash, open(hashfile, 'wb'))