# This uses the pearson networks to run the Boolean simulations.
# Each node is considered in a binary state, either on (1) or off (-1).
# The state of the node changes with each iteration.
# One of the nodes is used as an input signal, where it is on for 10 steps, then off for 10 steps.

import numpy as np
import pandas as pd
import time
import pickle
import networkx as nx
import json
import os
import re

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

def periodic_signal(len_period=20, iters=10):
    # iniciamos una señal periodica
    period = [1]*int(len_period/2) + [-1]*int(len_period/2) # dos listas concatenadas: una de 1 y otra de -1
    # 1 significa que está activada, -1 significa que está desactivada
    input_signal = period * int(iters / len(period)) + period[:np.mod(iters+1, len(period))]
    #np.mod(iters+1, len(period)) devuelve el resto de iters+1/len(period)
    #input_signal es entonces una lista (de iters+1 elementos) que repite el period
    #por ej si period=[1, 1, -1, -1], input_signal=[1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1]
    return input_signal

def add_noise(current_state, prob):
    partial_prob = np.asarray([prob] * len(current_state)) #array de 1d con len(current_state) elementos de valor prob
    flips = np.random.random(len(current_state)) # devuelve len(current_state) numeros aleatorios entre 0 y 1
    flips[flips <= partial_prob] = -1 #si es menor lo flipeamos (ponemos -1)
    flips[flips != -1] = 1 #el resto los dejamos igual (ponemos 1 porque multiplicar por 1 es quedarse igual)
    flips = flips.astype(int)
    new_state = current_state * flips #el nuevo estado es el anterior pero flipeando algunos elementos
    return new_state

node_dict = load_pickle('node_dictionary.pickle')
input_nodes = node_dict['GENETIC'] + node_dict['MOLECULAR'] + node_dict['PET'] + node_dict['MRI'] + node_dict['RISKFACTORS']
node_sort = input_nodes + node_dict['PHENOTYPE']

in_signal = True # decide if there is a manual signal driving the input
len_period = 20
noise_prob = 0.2 # 0, 0.05, 0.1, 0.15, 0.2 según el paper de Keith

# define the parameters
reps = 100 # number of times to run the simulation with new initial conditions for activation state
iters = 100 # number of time iterations to perform for the simulation

adj = load_network('../results_pearson/bignetwork_weighted_pearson.json') # upload la adj matrix de la pearson global network
npts = adj.shape[0]

if in_signal:
    input_signal = periodic_signal(len_period=len_period, iters=iters)
else:
    input_signal = np.zeros(iters+1)

since = time.time()
for i, in_node in enumerate(input_nodes):
    df_list = []
    for rep in range(reps):
        first_step = list(np.random.choice([-1, 1], size=npts)) # lista de -1 y 1 (random) y tamaño npts
        print(f'Node: {i+1:3d}/{len(input_nodes):3d}, Rep: {rep+1:3d}/{reps:3d}', end='\r', flush=True)
        active = pd.DataFrame(np.zeros((npts, iters+1)), index=node_sort, dtype=int)
        # dataframe donde las filas son cada nodo y las columnas son cada iteracion
        active.loc[:, 0] = first_step # la primera columna es el estado random del que partimos
        active.loc[in_node, :] = input_signal # la fila del nodo input es la señal que metemos, si no le metemos señal (input_signal es una
        #matriz de ceros), in_node será el único que tendrá un cero en la primera columna
        for step in range(iters): #en cada iteracion vamos llenando una columna del active
            current_state = active.loc[:, step]
            input_sum = adj.multiply(current_state, axis=0).sum(axis=0)
            # adj.multiply(current_state, axis=0) multiplica las filas de adj (los pearson coef) por current_state: fila 1 de adj * elemento 1 de current_state
            # como current_state es una lista de 1 y -1, el unico efecto es que pone negativos los que multiplica por -1
            # .sum(axis=0) suma todas las filas de una columna, input_sum tiene el nombre de las columnas y la suma
            next_step = input_sum.copy()
            # Asignamos -1 y 1 a cada nodo dependiendo de si la suma anterior es positiva o negativa
            # si es 0 lo dejamos como el paso anterior: next_step=current_state
            next_step[next_step > 0] = 1 # turn on nodes with positive incoming weight
            next_step[next_step < 0] = -1 # turn off nodes with negative incoming weight
            next_step[next_step == 0] = active.loc[list(next_step[next_step == 0].index), step] # keep previous step for 0 incoming weight
            next_step = next_step.astype(int)
            if noise_prob != 0:
                next_step = add_noise(next_step, noise_prob)
            active.loc[:, step+1] = next_step # update the next step
            if in_signal: #si la señal que metiamos eran 0 (nada), entonces el siguiente step es el que acabamos de calcular
                #si le metemos una señal nosotros, entonces el nodo que estamos perturbando lo cambiamos poniendole nuestra señal
                active.loc[in_node, step+1] = input_signal[step+1] # override input node with input signal
        #al final en cada repeticion tenemos un dataframe active de iters columnas
        df_list.append(active) #lista de los df active (numero de dataframes=reps)

    file = f'../results_pearson/network_paths/{in_node}/period_{len_period}/node_path_noise_{int(noise_prob * 100)}.pickle'
    create_folders(file)
    save_pickle(file, df_list) #guardamos el df_list para el input_node actual

print('')
time_elapsed = time.time() - since
print(f'Total time: {time_elapsed:.1f}s')