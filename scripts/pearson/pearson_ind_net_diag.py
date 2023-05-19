import numpy as np
import pandas as pd
import time
from itertools import combinations_with_replacement
import warnings
warnings.filterwarnings("ignore")
import pickle
import networkx as nx
import json

def save_pickle(file, thing, message=False):
    pickle_out = open(file, "wb")
    pickle.dump(thing, pickle_out)
    pickle_out.close()
    if message:
        print(f'new pickle created \'{file}\'')

def matrix_dim_from_triu(C, k=1):
    """
    k=1 assumes there are no self edges, change to k=0 if self edges are desired
    """
    dim = 0.5 * np.sqrt(1 + 8 * C)
    if k == 1:
        dim = int(0.5 + dim)
    elif k == 0:
        dim = int(-0.5 + dim)
    return dim

def adjacency_mat(coef, pval, df, alpha=0.05):
    """
    alpha: threshold for z-test to define an edge
    """
    nodes = list(df)
    ntups = len(coef)
    nfeats = matrix_dim_from_triu(ntups, k=0)
            
    sig_inds = np.where(pval < alpha)[0]
    tuples = list(combinations_with_replacement(range(nfeats), 2))
    sig_tups = [tuples[x] for x in sig_inds]
    
    adj_mat = np.zeros((nfeats, nfeats), dtype=int)
    weight_mat = np.zeros((nfeats, nfeats))
            
    for tupl in sig_tups:
        adj_mat[tupl] = 1
    adj_mat = np.tril(adj_mat.T, k=-1) + np.triu(adj_mat, k=1)
    weight_mat[np.triu_indices(nfeats, k=0)] = coef
    weight_mat.T[np.triu_indices(nfeats, k=0)] = coef
    weight_mat[np.isnan(weight_mat)] = 0
    weight_mat = adj_mat.astype(float) * weight_mat
    
    adj_df = pd.DataFrame(adj_mat, columns=nodes, index=nodes)
    weight_df = pd.DataFrame(weight_mat, columns=nodes, index=nodes)
    
    return [adj_df, weight_df]

def network_stats(G, df):
    """
    print some basic statistics of the network and the dataframe
    """
    degree = dict(G.degree)
    avg_deg = np.round(np.mean(list(degree.values())), decimals=1)
    avg_clst = nx.average_clustering(G) # average number of nodes
    print(f'Average Degree:     {avg_deg:>4.1f}')
    print(f'Number of Edges:    {G.number_of_edges():>4d}')
    print(f'Average Clustering: {avg_clst:>4.2f}')
    print('# Subjects:  ', df.shape[0])
    print('# Variables: ', df.shape[1])
    print('')

levels = ['GENETIC', 'MOLECULAR', 'PET', 'MRI', 'RISKFACTORS', 'PHENOTYPE']
diags = ['nl', 'mci', 'dementia']
print('Creating network')
since = time.time()
for level in levels:
    for diag in diags:
        print(level, diag)
        df = pd.read_csv(f'../datasets/{level}_{diag}.csv')
        df.set_index('RID', inplace=True)
        coef = np.load(f'../results_pearson/{level}_{diag}_pearson_array.npy')
        pval = np.load(f'../results_pearson/{level}_{diag}_pval_array.npy')
        mat = adjacency_mat(coef, pval, df)
        for i, name in enumerate(["unweighted", "weighted"]):
            G = nx.from_pandas_adjacency(mat[i]) # create networkx object from adjacency matrix
            G_json = nx.readwrite.json_graph.cytoscape.cytoscape_data(G) # create json (dictionary object) to store network info
            adj_file = f'../results_pearson/{level}_{diag}_adj_df_{name}_pearson.pickle'
            net_file = f'../results_pearson/{level}_{diag}_network_{name}_pearson.json'
            save_pickle(adj_file, mat[i]) # save the adjacency matrix (pandas dataframe)
            json.dump(G_json, open(net_file, 'w')) # file to manipulate later if needed (change names, add attributes, etc.)
        print(adj_file)
        print(net_file)
        network_stats(G, df)
time_elapsed = time.time() - since
print(f'Total time: {time_elapsed:.1f}s')