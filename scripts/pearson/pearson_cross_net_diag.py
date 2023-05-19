import numpy as np
import pandas as pd
import time
from itertools import combinations_with_replacement, combinations
import warnings
warnings.filterwarnings("ignore")

levels = ['GENETIC', 'MOLECULAR', 'PET', 'MRI', 'RISKFACTORS', 'PHENOTYPE']
diags = ['nl', 'mci', 'dementia']
combos = list(combinations(levels, 2))

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

print('Creating cross networks')
since = time.time()
for level1, level2 in combos:
    for diag in diags:
        print(level1, level2, diag)
        coef = np.load(f'../results_pearson/pearson_array_{level1}_{level2}_{diag}.npy')
        pval = np.load(f'../results_pearson/pval_array_{level1}_{level2}_{diag}.npy')
        df1 = pd.read_csv(f'../datasets/{level1}_{diag}.csv')
        df1.set_index('RID', inplace=True)
        df2 = pd.read_csv(f'../datasets/{level2}_{diag}.csv')
        df2.set_index('RID', inplace=True)
        nodes1, nodes2 = list(df1), list(df2)

        df_joined = join_df(df1, df2, level1, level2)
                
        mat = adjacency_mat(coef, pval, df_joined)
        
        adj1 = load_pickle(f'../results_pearson/{level1}_{diag}_adj_df_weighted_pearson.pickle') #las networks calculadas antes
        adj2 = load_pickle(f'../results_pearson/{level2}_{diag}_adj_df_weighted_pearson.pickle')

        # para incluir los edges dentro de las networks de cada grupo (osea para tener la figura final tocha)
        # use edges from individual level networks (higher sample sizes)
        individual_edges(nodes1, adj1, mat[0], mat[1])
        individual_edges(nodes2, adj2, mat[0], mat[1])

        for i, name in enumerate(["unweighted", "weighted"]):
            G = nx.from_pandas_adjacency(mat[i])
            G_json = nx.readwrite.json_graph.cytoscape.cytoscape_data(G)
            adj_file = f'../results_pearson/{level1}_{level2}_{diag}_adj_df_{name}_pearson.pickle'
            net_file = f'../results_pearson/{level1}_{level2}_{diag}_network_{name}_pearson.json'
            save_pickle(adj_file, mat[i])
            json.dump(G_json, open(net_file, 'w'))
                
        print(adj_file) #printemos solo las cosas del weighted
        print(net_file)
        network_stats(G, df_joined)
time_elapsed = time.time() - since
print(f'Total time: {time_elapsed:.1f}s')