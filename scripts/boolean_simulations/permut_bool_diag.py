import networkx as nx
from networkx.algorithms import isomorphism
import random
import pickle
import json
import time
from datetime import timedelta
import collections
import itertools
import numpy as np
import pandas as pd
import scipy.signal
from joblib import Parallel, delayed

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

class EdgeSwapGraph(nx.Graph):
# Copyright (c) 2011-2012 Christopher D. Lasher
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
  """An interaction graph which can produce a "random" graph from
  itself by an iterative edge-swap technique that preserves the degree
  distribution of the original graph.
  """
  def randomize_by_edge_swaps(self, num_iterations):
    """Randomizes the graph by swapping edges in such a way that
    preserves the degree distribution of the original graph.
    The underlying idea stems from the following. Say we have this
    original formation of edges in the original graph:
        head1   head2
          |       |
          |       |
          |       |
        tail1   tail2
    Then we wish to swap the edges between these four nodes as one
    of the two following possibilities:
        head1   head2       head1---head2
              \ /
              X
              / \
        tail1   tail2       tail1---tail2
    We approach both by following through the first of the two
    possibilities, but before committing to the edge creation, give
    a chance that we flip the nodes `head1` and `tail1`.
    See the following references for the algorithm:
    - F. Viger and M. Latapy, "Efficient and Simple Generation of
      Random Simple Connected Graphs with Prescribed Degree
      Sequence," Computing and Combinatorics, 2005.
    - M. Mihail and N.K. Vishnoi, "On Generating Graphs with
      Prescribed Vertex Degrees for Complex Network Modeling,"
      ARACNE 2002, 2002, pp. 1–11.
    - R. Milo, N. Kashtan, S. Itzkovitz, M.E.J. Newman, and U.
      Alon, "On the uniform generation of random graphs with
      prescribed degree sequences," cond-mat/0312028, Dec. 2003.
    :Parameters:
    - `num_iterations`: the number of iterations for edge swapping
      to perform; this value will be multiplied by the number of
      edges in the graph to get the total number of iterations
    """
    newgraph = self.copy()
    edge_list = list(newgraph.edges()) # KK change, added list()
    num_edges = len(edge_list)
    total_iterations = num_edges * num_iterations

    for _ in range(total_iterations):
      rand_index1 = int(round(random.random() * (num_edges - 1)))
      rand_index2 = int(round(random.random() * (num_edges - 1)))
      original_edge1 = edge_list[rand_index1]
      original_edge2 = edge_list[rand_index2]
      head1, tail1 = original_edge1
      head2, tail2 = original_edge2

      # Flip a coin to see if we should swap head1 and tail1 for
      # the connections
      if random.random() >= 0.5:
        head1, tail1 = tail1, head1

      # The plan now is to pair head1 with tail2, and head2 with
      # tail1
      #
      # To avoid self-loops in the graph, we have to check that,
      # by pairing head1 with tail2 (respectively, head2 with
      # tail1) that head1 and tail2 are not actually the same
      # node. For example, suppose we have the edges (a, b) and
      # (b, c) to swap.
      #
      #   b
      #  / \
      # a   c
      #
      # We would have new edges (a, c) and (b, b) if we didn't do
      # this check.

      if head1 == tail2 or head2 == tail1:
        continue

      # Trying to avoid multiple edges between same pair of nodes;
      # for example, suppose we had the following
      #
      # a   c
      # |*  |           | original edge pair being looked at
      # | * |
      # |  *|           * existing edge, not in the pair
      # b   d
      #
      # Then we might accidentally create yet another (a, d) edge.
      # Note that this also solves the case of the following,
      # missed by the first check, for edges (a, b) and (a, c)
      #
      #   a
      #  / \
      # b   c
      #
      # These edges already exist.

      if newgraph.has_edge(head1, tail2) or newgraph.has_edge(head2, tail1):
        continue

      # Suceeded checks, perform the swap
      original_edge1_data = newgraph[head1][tail1]
      original_edge2_data = newgraph[head2][tail2]

      newgraph.remove_edges_from((original_edge1, original_edge2))

      new_edge1 = (head1, tail2, original_edge1_data)
      new_edge2 = (head2, tail1, original_edge2_data)
      newgraph.add_edges_from((new_edge1, new_edge2))

      # Now update the entries at the indices randomly selected
      edge_list[rand_index1] = (head1, tail2)
      edge_list[rand_index2] = (head2, tail1)

    assert len(newgraph.edges()) == num_edges
    return newgraph

def PermutNetwork(G, levels, combos, level_nodes):
  # FOR BETWEEN LAYERS --> Lo tratamos como una bipartite network: solo tenemos en cuenta edges entre layers
  for i, (level1, level2) in enumerate(combos):
    # print(f'Combo:{i+1:2d}/{len(combos)}', end='\r', flush=True)
    C = G.subgraph(level_nodes[level1] + level_nodes[level2]).copy()
    edgesOld = [x for x in C.edges() if x[0] in level_nodes[level1] and x[1] in level_nodes[level2]]
    edgesOld = edgesOld + [x[::-1] for x in C.edges() if x[1] in level_nodes[level1] and x[0] in level_nodes[level2]]
    sources = [x[0] for x in edgesOld]
    targets = [x[1] for x in edgesOld]
    weights = [C[x[0]][x[1]] for x in edgesOld]

    # Permutamos los edges * --> en una bipartite network es tan sencillo como desordenar los nodos de la segunda red
    targetRand = random.sample(targets, len(targets)) # reordenamos los targets aleatoriamente
    edgesNew = list(zip(sources, targetRand, weights)) # nuevos edges entre los source y los targets reordenados
    sub_edges = [(x, y) for x, y, _ in edgesNew]
    counts = collections.Counter(sub_edges)
    maxCount = max(counts.values())

    while maxCount > 1:
      singleEdges = [(edgesNew[sub_edges.index(key)]) for key, val in counts.items() if val == 1]
      multiEdges = [(edgesNew[sub_edges.index(key)]) for key, val in counts.items() if val > 1] # se repite el edge (el mismo par de nodos tiene mas de un edge)
      for multEdge in multiEdges: # si aparecen multiedges en la bipartite network
        singleEdges2 = singleEdges.copy()
        redundant = [0]
        while len(redundant) > 0 and len(singleEdges2) > 0:
          singEdge = random.sample(singleEdges2, 1)[0] # coge aleatoriamente un edge
          newMult = (multEdge[0], singEdge[1], multEdge[2])
          newSing = (singEdge[0], multEdge[1], singEdge[2])
          list_sub = [(x,y) for x, y, _ in [newMult, newSing]]
          redundant = [x for x in sub_edges if x in list_sub]
          singleEdges2.remove(singEdge)
          # print(f'{len(singleEdges2):6d}/{len(singleEdges)}', end='\r')
        if len(singleEdges2) <= 0:
          break # break el for loop y pasamos al siguiente if **
        else:
          singleEdges[singleEdges.index(singEdge)] = newSing
          edgesNew[edgesNew.index(singEdge)] = newSing
          edgesNew[edgesNew.index(multEdge)] = newMult
      
      # **
      if len(singleEdges2) <= 0:
        # hacemos otra permutacion para la bipartite network
        # permutamos los edges como en * y volvemos al while como antes
        targetRand = random.sample(targets, len(targets))
        edgesNew = list(zip(sources, targetRand, weights))
        sub_edges = [(x, y) for x, y, _ in edgesNew]
        counts = collections.Counter(sub_edges)
        maxCount = max(counts.values())
      else:
        sub_edges = [(x, y) for x, y, _ in edgesNew]
        counts = collections.Counter(sub_edges)
        maxCount = max(counts.values())

    G.remove_edges_from(edgesOld) #lo he añadido yo
    G.add_edges_from(edgesNew) #lo he añadido yo


  # FOR SINGLE LAYERS
  # Codigo de Christopher D. Lasher
  for i, level in enumerate(levels):
    # print(f'Layer: {i+1:1d}/{len(levels)}', end='\r', flush=True)
    S = G.subgraph(level_nodes[level]).copy()
    N = EdgeSwapGraph(S).randomize_by_edge_swaps(num_iterations=1000)
    edgesOld = list(S.edges(data=True))
    edgesNew = list(N.edges(data=True))
    G.remove_edges_from(edgesOld)
    G.add_edges_from(edgesNew)

  return G

def periodic_signal(len_period=20, iters=10):
    period = [1]*int(len_period/2) + [-1]*int(len_period/2)
    input_signal = period * int(iters / len(period)) + period[:np.mod(iters+1, len(period))]
    return input_signal

def add_noise(current_state, prob):
    partial_prob = np.asarray([prob] * len(current_state))
    flips = np.random.random(len(current_state))
    flips[flips <= partial_prob] = -1
    flips[flips != -1] = 1
    flips = flips.astype(int)
    new_state = current_state * flips
    return new_state

def k_shortest_paths(G, source, target, k, weight=None):
    return list(itertools.islice(nx.algorithms.simple_paths.shortest_simple_paths(G, source, target, weight=weight), k))

def BooleanSimulation(G, input_nodes, node_sort, targets, len_period=20, nreps=3, iters=3, noise_prob=0.05, in_signal=True, n_paths=10):
  """Performs the boolean simulations, calculates cross correlation coefficient and the top 10 shortest paths.
  """
  adj = nx.convert_matrix.to_pandas_adjacency(G, nodelist=list(G.nodes()))
  npts = adj.shape[0]
  tuples = list(itertools.combinations_with_replacement(range(len(node_sort)), 2))
  
  if in_signal:
    input_signal = periodic_signal(len_period=len_period, iters=iters)
  else:
    input_signal = np.zeros(iters+1)
  
  path_dict = {}
  # Boolean simulation
  print('Start boolean simulation')
  for i, in_node in enumerate(input_nodes):
    cross = np.ones((nreps, npts, npts))
    for rep in range(nreps):
      # print(f'Node: {i+1:3d}/{len(input_nodes)}, Rep: {rep+1:3d}/{nreps}', end='\r', flush=True)
      first_step = list(np.random.choice([-1, 1], size=npts))
      active = pd.DataFrame(np.zeros((npts, iters+1)), index=node_sort, dtype=int)
      active.loc[:, 0] = first_step
      active.loc[in_node, :] = input_signal
      for step in range(iters):
        current_state = active.loc[:, step]
        input_sum = adj.multiply(current_state, axis=0).sum(axis=0)
        next_step = input_sum.copy()
        next_step[next_step > 0] = 1
        next_step[next_step < 0] = -1
        next_step[next_step == 0] = active.loc[list(next_step[next_step == 0].index), step]
        next_step = next_step.astype(int)
        if noise_prob != 0:
          next_step = add_noise(next_step, noise_prob)
        active.loc[:, step+1] = next_step
        if in_signal:
          active.loc[in_node, step+1] = input_signal[step+1]
      # Cross correlation coefficient
      rep_path = active.values
      for tupl in tuples:
        x = rep_path[tupl[0]]
        y = rep_path[tupl[1]]
        cross[rep][tupl] = np.max(np.abs(scipy.signal.correlate(x, y)))
        cross[rep][tupl[::-1]] = np.max(np.abs(scipy.signal.correlate(x, y)))
      cross[rep] /= np.max(cross[rep])
    cross_mean = pd.DataFrame(np.mean(cross, axis=0), index=node_sort, columns=node_sort)
    print('Finished boolean simulation')

    # Shortest path
    print('Start shortest path')
    adj_cross = adj.copy()
    adj_cross[adj_cross!=0] = 1
    adj_cross = adj_cross.loc[node_sort, node_sort]
    adj_cross = adj_cross * cross_mean
    adj_cross[adj_cross>0] = 1 / adj_cross[adj_cross>0]
    net_cross = nx.convert_matrix.from_pandas_adjacency(adj_cross)
    
    path_dict[in_node] = {}
    for target in targets:
      try:
        path_dict[in_node][target] = k_shortest_paths(net_cross, in_node, target, n_paths, weight='weight')
      except:
        path_dict[in_node][target] = []
  
  return path_dict

def are_graphs_isomorphic(graph1, graph2):
    # Pruning inicial
    if graph1.number_of_nodes() != graph2.number_of_nodes() or \
            graph1.number_of_edges() != graph2.number_of_edges():
        return False

    # Ordenar nodos por grado
    nodes1 = sorted(graph1.nodes(), key=lambda x: graph1.degree(x))
    nodes2 = sorted(graph2.nodes(), key=lambda x: graph2.degree(x))
    mapping = {nodes1[i]: nodes2[i] for i in range(len(nodes1))}
    graph1_reorder = nx.relabel_nodes(graph1, mapping)

    vf2 = isomorphism.vf2userfunc.GraphMatcher(graph1_reorder, graph2)
    return vf2.is_isomorphic()

def parallel_permutation(permut, diag, num, G, levels, combos, node_dict):
  input_nodes = node_dict['GENETIC'] + node_dict['MOLECULAR'] + node_dict['PET'] + node_dict['MRI'] + node_dict['RISKFACTORS']
  targets = node_dict['PHENOTYPE']
  node_sort = input_nodes + targets
  print('Upload network')
  P = G.copy()
  print('Start permuting')
  step = 0
  while are_graphs_isomorphic(G, P): #si las redes son iguales se vuelve a permutar
    print(f'{step:5d}', end='\r')
    P = PermutNetwork(P, levels, combos, node_dict)
    step += 1
  print('Finished permutation')
  P_json = nx.readwrite.json_graph.cytoscape.cytoscape_data(P)
  json.dump(P_json, open(f'../results_pearson/negative_controls/network_{diag}_{num}_permuted{permut}.json', 'w'))
  print('Network saved!')
  paths = BooleanSimulation(P, input_nodes, node_sort, targets)
  save_pickle(f'../results_pearson/negative_controls/top_pathways_{diag}_{num}_permuted{permut}.pickle', paths)
  print('Top paths saved!')

levels = ['GENETIC', 'MOLECULAR', 'PET', 'MRI', 'RISKFACTORS', 'PHENOTYPE']
diags = ['nl', 'mci', 'dementia']
combos = list(itertools.combinations(levels, 2))
num_permutations = 100

for num in range(4):
  since = time.time()
  for diag in diags:
    node_dict = load_pickle(f'../results_pearson/para_controls/node_dictionary_{diag}_{num}.pickle')
    G = load_networkx(f'../results_pearson/para_controls/network_{diag}_{num}.json')
    for u, v in G.edges():
      del G[u][v]['source']
      del G[u][v]['target']

    # Parallel(n_jobs=num_permutations)(delayed(parallel_permutation)(j) for j in range(num_permutations))
    for j in range(num_permutations):
      parallel_permutation(j, diag, num, G, levels, combos, node_dict)
  time_elapsed = time.time() - since
  print(f'Total time: {str(timedelta(seconds=time_elapsed))}')