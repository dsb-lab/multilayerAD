import pickle
import json
import itertools
import networkx as nx

def load_pickle(file):
    """
    load a pickle file into a local variable
    can be a pandas dataframe or dictionary
    """
    pickle_in = open(file, "rb")
    pickle_file = pickle.load(pickle_in)
    pickle_in.close()
    return pickle_file

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

def load_network(file, node_list=None):
    G = load_networkx(file)
    if node_list == None:
        node_list = list(G.nodes())
    adj = nx.convert_matrix.to_pandas_adjacency(G, nodelist=node_list) # gives the graph adjacency matrix as a pandas df
    return adj

# Creating networks from cross correlation files
def create_cross_corr_net(diag, source, len_period, noise_prob, node_list=None):
    adj = load_network(f'../results_pearson/bignetwork_{diag}_weighted_pearson.json', node_list=node_list)
    cross_dict = load_pickle(f'../results_pearson/network_paths/{source}/period_{len_period}/cross_corr_{diag}_noise_{int(noise_prob*100)}.pickle') #guardabamos mean y std
    cross_mean = cross_dict['mean']
    adj[adj!=0] = 1
    adj = adj.loc[node_list, node_list] # reordenamos para la multiplicacion de la siguiente linea
    adj_cross = adj * cross_mean
    # todo lo anterior es para sustituir el valor de la adj matrix por el del cross-corr
    # para los elementos de la matriz que son distintos de 0:
    adj_cross[adj_cross>0] = 1 / adj_cross[adj_cross>0] # shortest path algorithm uses smaller weights, want higher cross correlations to have smaller weights
    net_cross = nx.convert_matrix.from_pandas_adjacency(adj_cross)
    return net_cross

def k_shortest_paths(G, source, target, k, weight=None):
    #se hace el shortest path algorithm k veces
    return list(itertools.islice(nx.algorithms.simple_paths.shortest_simple_paths(G, source, target, weight=weight), k))

node_dict = load_pickle('node_dictionary.pickle')
sources = node_dict['GENETIC'] + node_dict['MOLECULAR'] + node_dict['PET'] + node_dict['MRI'] + node_dict['RISKFACTORS']
targets = node_dict['PHENOTYPE']
node_sort = sources + targets
diags = ['nl', 'mci', 'dementia']

len_period = 20
noise_prob = 0.05 #escoger el que dé valores de cross correlation "intermedios"-> que no sean casi todos 1
n_paths = 10

path_dict = {}
for h, diag in enumerate(diags):
    path_dict[diag] = {}
    for i, source in enumerate(sources):
        G = create_cross_corr_net(diag, source, len_period, noise_prob, node_list=node_sort)
        path_dict[diag][source] = {}
        for j, target in enumerate(targets):
            print(f'{h+1:1d}/{len(diags):1d} {diag:>8s}, {i+1:2d}/{len(sources):2d} {source:<20}  --> {target:>20}', end='\r')
            try:
                path_dict[diag][source][target] = k_shortest_paths(G, source, target, n_paths, weight='weight')
            except:
                path_dict[diag][source][target] = []
save_pickle('../results_pearson/network_paths/top_pathways_diag.pickle', path_dict)