import numpy as np
import pandas as pd
import time
from itertools import combinations_with_replacement, combinations
import warnings
warnings.filterwarnings("ignore")

levels = ['GENETIC', 'MOLECULAR', 'PET', 'MRI', 'RISKFACTORS', 'PHENOTYPE']
combos = list(combinations(levels, 2))

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

def individual_edges(nodes, adj_indv_df, adj_df, weight_df):
    # adj_indv_df es la matriz con pesos
    pairs = list(combinations(nodes, 2))
    for pair in pairs:
        adj_df.loc[pair] = np.ceil(adj_indv_df.loc[pair]).astype(int)
        # np.ceil devuelve una matriz de los ceils de cada elemento
        # The ceil of the scalar x is the smallest integer i, such that i >= x. Ej: 1.7 -> 2, 0.2 -> 1
        # En nuestro caso, los elementos de adj1 son valores entre 0-1. Los 0 se quedan igual y los > 0 se hacen 1
        # Queremos esto porque cualquier correlacion es válida, ya que han pasado el permutation test
        weight_df.loc[pair] = adj_indv_df.loc[pair]

def load_pickle(file):
    """
    load a pickle file into a local variable
    can be a pandas dataframe or dictionary
    """
    pickle_in = open(file, "rb")
    pickle_file = pickle.load(pickle_in)
    pickle_in.close()
    return pickle_file

def join_df(df1, df2, level1, level2):
    """
    find subjects that overlap between the two dataframes
    create new dataframe that combines the two level dataframes
    """
    inds = [x for x in df1.index if x in df2.index]
    print(f'# points for {level1}: {df1.shape[0]}')
    print(f'# points for {level2}: {df2.shape[0]}')
    print(f'# points overlap: {len(inds)}')
    df1_cut = df1.loc[inds]
    df2_cut = df2.loc[inds]
    df_joined = pd.concat([df1_cut, df2_cut], axis=1)
    
    return df_joined

for level1, level2 in combos:
    print(level1, level2)
    coef_real = np.load(f'../results_mi/mi_array_{level1}_{level2}.npy')
    coef_perm = np.load(f'../results_mi/mi_perms_{level1}_{level2}.npy')
    df1 = pd.read_csv(f'../datasets/{level1}.csv')
    df1.set_index('RID', inplace=True)
    df2 = pd.read_csv(f'../datasets/{level2}.csv')
    df2.set_index('RID', inplace=True)
    nodes1, nodes2 = list(df1), list(df2)

    df_joined = join_df(df1, df2, level1, level2)
            
    mat = adjacency_mat(coef_real, coef_perm, df_joined)
    
    adj1 = load_pickle(f'../results_mi/{level1}_adj_df_weighted_mi.pickle') #las networks calculadas antes
    adj2 = load_pickle(f'../results_mi/{level2}_adj_df_weighted_mi.pickle')

    # para incluir los edges dentro de las networks de cada grupo (osea para tener la figura final tocha)
    # use edges from individual level networks (higher sample sizes)
    individual_edges(nodes1, adj1, mat[0], mat[1])
    individual_edges(nodes2, adj2, mat[0], mat[1])

    for i, name in enumerate(["unweighted", "weighted"]):
        G = nx.from_pandas_adjacency(mat[i])
        G_json = nx.readwrite.json_graph.cytoscape.cytoscape_data(G)
        adj_file = f'../results_mi/{level1}_{level2}_adj_df_{name}_mi.pickle'
        net_file = f'../results_mi/{level1}_{level2}_network_{name}_mi.json'
        save_pickle(adj_file, mat[i])
        json.dump(G_json, open(net_file, 'w'))
            
    print(adj_file)
    print(net_file)
    network_stats(G, df_joined)