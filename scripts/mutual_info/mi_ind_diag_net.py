import numpy as np
import pandas as pd
import time
from itertools import combinations_with_replacement
import warnings
warnings.filterwarnings("ignore")

levels = ['GENETIC', 'MOLECULAR', 'PET', 'MRI', 'RISKFACTORS', 'PHENOTYPE']
add_degree = [0, 1000, 2000, 3000, 5000, 4000]
diags = ['nl', 'mci', 'dementia']

import pickle
from scipy.stats import norm
from sklearn.metrics import normalized_mutual_info_score
import networkx as nx
import json

def save_pickle(file, thing, message=False):
    pickle_out = open(file, "wb")
    pickle.dump(thing, pickle_out)
    pickle_out.close()
    if message:
        print(f'new pickle created \'{file}\'')

def coef_norm_test(coef_real, coef_perm):
    """
    compares the real MI value to the distribution of MI scores from the permutations
    returns the p-value and the z-score
    """
    pmean = np.nanmean(coef_perm) # hace el mean pero ignorando nans
    pstd = np.nanstd(coef_perm) # hace std pero ignorando nans
    if pstd == 0:
        # in case the permutations produced the same MI score for every repetition
        pval = np.nan
    else:
        zscore = (coef_real - pmean) / pstd
        pval = 1 - norm.cdf(abs(zscore))
    # A normal cumulative distribution function (CDF) return the percentage of the normal distribution function that is less
    # than or equal to the random variable specified.
    return pval

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

def adjacency_mat(coef_real, coef_perm, df, alpha=0.05):
    """
    alpha: threshold for z-test to define an edge
    """
    nodes = list(df)
    ntups = coef_perm.shape[1] # shape = number of repetitions, number of pairs
    nfeats = matrix_dim_from_triu(ntups, k=0) # calculate original number of variables
    
    pval = np.zeros(ntups) # define empty matrix to store p-values for z-testing

    for i in range(ntups):
        pval[i] = coef_norm_test(coef_real[i], coef_perm[:, i])
    pval[np.isnan(pval)] = 1 # in case there are errors, p-value of 1 is insignificant
            
    sig_inds = np.where(pval < alpha)[0] # find indices variable pairs with significant p-values
    tuples = list(combinations_with_replacement(range(nfeats), 2)) # create list of tuples (variable pairs)
    sig_tups = [tuples[x] for x in sig_inds] # list of actual significant tuples (variable pairs)
    
    adj_mat = np.zeros((nfeats, nfeats), dtype=int) # empty adjacency matrix for binary network (1: edge, 0: no edge)
    # guardamos la adj_mat por si la necesitamos en algún momento pero no es tremendamente util
    weight_mat = np.zeros((nfeats, nfeats)) # empty adj mat for weighted network (MI score for edge: 0-1, 0: no edge)
            
    for tupl in sig_tups:
        if tupl[0] != tupl[1]:
            adj_mat[tupl] = 1
            z = np.array(df.iloc[:, list(tupl)], dtype=float).T
            nan_inds = np.where(np.isnan(z))[1]
            z = np.delete(z, nan_inds, axis=1)
            weight_mat[tupl] = normalized_mutual_info_score(*z) # use normalized score between 0-1
        
    # fill in the lower triangle of the two adjacency matrices
    adj_mat = np.triu(adj_mat, k=1).T + np.triu(adj_mat, k=1)
    # en el triangulo sup tiene los valores MI (la diagonal son 0). k=1 en triu_indices hace que no te coja la diagonal
    # la matriz que se obtiene es simetrica
    weight_mat = np.triu(weight_mat, k=1).T + np.triu(weight_mat, k=1)
        
    # create a pandas dataframe from the unweighted adjacency matrix, this way we can add the names of the nodes
    adj_df = pd.DataFrame(adj_mat, columns=nodes, index=nodes)
    # create a pandas dataframe from the weighted adjacency matrix
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

print('Creating network')
since = time.time()
for deg, level in enumerate(levels):
    for diag in diags:
        print(level, diag)
        df = pd.read_csv(f'../datasets/{level}_{diag}.csv')
        df.set_index('RID', inplace=True)
        coef_real = np.load(f'../results_mi/{level}_{diag}_mi_array.npy')
        coef_perm = np.load(f'../results_mi/{level}_{diag}_mi_perms.npy')
        mat = adjacency_mat(coef_real, coef_perm, df)
        for i, name in enumerate(["unweighted", "weighted"]):
            G = nx.from_pandas_adjacency(mat[i]) # create networkx object from adjacency matrix
            values = dict(G.degree)
            values = {k: v+add_degree[deg] for k,v in values.items()}
            nx.set_node_attributes(G, values, name='node_degree')
            G_json = nx.readwrite.json_graph.cytoscape.cytoscape_data(G) # create json (dictionary object) to store network info
            adj_file = f'../results_mi/{level}_{diag}_adj_df_{name}_mi.pickle'
            net_file = f'../results_mi/{level}_{diag}_network_{name}_mi.json'
            save_pickle(adj_file, mat[i]) # save the adjacency matrix (pandas dataframe)
            json.dump(G_json, open(net_file, 'w')) # file to manipulate later if needed (change names, add attributes, etc.)
        print(adj_file)
        print(net_file)
        network_stats(G, df)
time_elapsed = time.time() - since
print(f'Total time: {time_elapsed:.1f}s')