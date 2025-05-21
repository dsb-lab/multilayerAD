import numpy as np
import pandas as pd
import pickle
import scipy.signal
import os
import re
import time
from itertools import combinations_with_replacement

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

def create_folders(file):
    slash_inds = [match.start() for match in re.finditer('/', file)]
    for ind in slash_inds:
        sub = file[:ind+1]
        if not os.path.isdir(sub):
            os.mkdir(sub)
    
node_dict = load_pickle('node_dictionary.pickle')
input_nodes = node_dict['GENETIC'] + node_dict['MOLECULAR'] + node_dict['PET'] + node_dict['MRI'] + node_dict['RISKFACTORS']
node_sort = input_nodes + node_dict['PHENOTYPE']

len_period = 20
noise_prob = 0 # 0, 0.05, 0.1, 0.15, 0.2

since = time.time()
for k, in_node in enumerate(input_nodes):
    print(f'Input: {k+1:02d}/{len(input_nodes):02d}', end='\r')
    path_list = load_pickle(f'../results_pearson/network_paths/{in_node}/period_{len_period}/node_path_noise_{int(noise_prob*100)}.pickle')
    # recordamos que el path_list es una lista para un determinado nodo input en la que cada elemento es una rep

    # define the pairs to calculate the cross correlation between
    tuples = list(combinations_with_replacement(range(len(node_sort)), 2))
    nreps = len(path_list)
    nfeats, niters = path_list[0].shape
    cross = np.ones((nreps, nfeats, nfeats))

    # actually calculate the cross-corr
    # el cross-corr se calcula entre dos nodos y nos da una idea de como de relacionados estan (si la variacion de uno afecta al otro)
    for rep in range(nreps): #cogemos para cada repeticion
        rep_path = path_list[rep].values
        #coge solo los valores, por lo que nos queda una lista en la que cada elemento es una lista que contiene la fila
        for tupl in tuples:
            x = rep_path[tupl[0]] #fila del nodo 1
            y = rep_path[tupl[1]] #fila del nodo 2
            cross[rep][tupl] = np.max(np.abs(scipy.signal.correlate(x, y)))
            #usamos scipy.signal.correlate(x, y) y no np.correlate(x,y,mode='full') porque la de scipy usa fft que es mejor para arrays grandes
            #el cross-corr es un valor entre -1 (anticorr) y 1 (corr), pero a nosotros solo nos interesa si estan correlacionadas o no, hacemos np.abs
            #del paper de keith: the maximum value of the temporal cross correlation across all lag times -> PREGUNTAR
            # con np.max cogemos el valor más grande de todos los calculados
            cross[rep][tupl[::-1]] = np.max(np.abs(scipy.signal.correlate(x, y))) # porque el coef es igual sea la tupla (1,2) que (2,1)
        cross[rep] /= np.max(cross[rep]) # normalize the cross correlation values between 0 and 1
        
    cross_mean = pd.DataFrame(np.mean(cross, axis=0), index=node_sort, columns=node_sort)
    cross_std = pd.DataFrame(np.std(cross, axis=0), index=node_sort, columns=node_sort)
    cross_dict = {'mean':cross_mean, 'std':cross_std}

    # define file, create directories if needed, save results
    file = f'../results_pearson/network_paths/{in_node}/period_{len_period}/cross_corr_noise_{int(noise_prob*100)}.pickle'
    create_folders(file)
    save_pickle(file, cross_dict)
time_elapsed = time.time() - since
print(f'Total time: {time_elapsed:.1f}s')