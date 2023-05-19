import pickle
import networkx as nx
import json

def save_pickle(file, thing, message=False):
    pickle_out = open(file, "wb")
    pickle.dump(thing, pickle_out)
    pickle_out.close()
    if message:
        print(f'new pickle created \'{file}\'')

def load_networkx(file):
    with open(file) as infile:
        network = json.load(infile)
    G = nx.readwrite.json_graph.cytoscape_graph(network)
    return G

# load networks
G_gene = load_networkx('../results_mi/GENETIC_network_weighted_mi.json')
G_phen = load_networkx('../results_mi/PHENOTYPE_network_weighted_mi.json')
Gtot = load_networkx('../results_mi/multilevelnetwork_weighted_mi.json')
for _, _, d in Gtot.edges(data=True):
    d['weight'] = 1/d['weight']

sources = list(sorted(G_gene.nodes()))
targets = list(sorted(G_phen.nodes()))

path_dict = {}
for i, source in enumerate(sources):
    path_dict[source] = {}
    for j, target in enumerate(targets):
        print(f'{i+1:02d}/{len(sources):02d} {source} --> {target}', end='\r')
        try:
            path_dict[source][target] = nx.shortest_path(Gtot, source, target, weight='weight')
        except: # if there is no path connecting the two nodes
            path_dict[source][target] = []

save_pickle('../results_mi/top_pathways_multinet_mi.pickle', path_dict)