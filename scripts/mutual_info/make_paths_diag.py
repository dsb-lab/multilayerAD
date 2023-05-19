import pickle
import networkx as nx
import itertools
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

diags = ['nl', 'mci', 'dementia']
path_dict = {}
for h, diag in enumerate(diags):
    # load networks
    G_gene = load_networkx(f'../results_mi/GENETIC_{diag}_network_weighted_mi.json')
    G_phen = load_networkx(f'../results_mi/PHENOTYPE_{diag}_network_weighted_mi.json')
    G_gene_biom = load_networkx(f'../results_mi/GENETIC_MOLECULAR_{diag}_network_weighted_mi.json')
    G_gene_risk = load_networkx(f'../results_mi/GENETIC_RISKFACTORS_{diag}_network_weighted_mi.json')
    G_biom_pet = load_networkx(f'../results_mi/MOLECULAR_PET_{diag}_network_weighted_mi.json')
    G_biom_mri = load_networkx(f'../results_mi/MOLECULAR_MRI_{diag}_network_weighted_mi.json')
    G_biom_risk = load_networkx(f'../results_mi/MOLECULAR_RISKFACTORS_{diag}_network_weighted_mi.json')
    G_mri_pet = load_networkx(f'../results_mi/PET_MRI_{diag}_network_weighted_mi.json')
    G_pet_phen = load_networkx(f'../results_mi/PET_PHENOTYPE_{diag}_network_weighted_mi.json')
    G_mri_phen = load_networkx(f'../results_mi/MRI_PHENOTYPE_{diag}_network_weighted_mi.json')
    G_risk_pet = load_networkx(f'../results_mi/PET_RISKFACTORS_{diag}_network_weighted_mi.json')
    G_risk_mri = load_networkx(f'../results_mi/MRI_RISKFACTORS_{diag}_network_weighted_mi.json')
    G_risk_phen = load_networkx(f'../results_mi/RISKFACTORS_PHENOTYPE_{diag}_network_weighted_mi.json')
    # Multilevel network:
    Gtot = nx.compose(nx.compose(nx.compose(nx.compose(nx.compose(nx.compose(nx.compose(nx.compose(nx.compose(G_gene_biom, G_biom_pet), G_biom_mri), G_mri_pet), G_mri_phen), G_pet_phen), G_risk_phen), G_biom_risk), G_risk_pet), G_risk_mri)
    for _, _, d in Gtot.edges(data=True):
        d['weight'] = 1/d['weight']

    sources = list(sorted(G_gene.nodes()))
    targets = list(sorted(G_phen.nodes()))

    path_dict[diag] = {}
    for i, source in enumerate(sources):
        path_dict[diag][source] = {}
        for j, target in enumerate(targets):
            print(f'{h+1:d}/{len(diags):d} {diag}, {i+1:02d}/{len(sources):02d} {source} --> {target}', end='\r')
            try:
                path_dict[diag][source][target] = nx.shortest_path(Gtot, source, target, weight='weight')
            except: # if there is no path connecting the two nodes
                path_dict[diag][source][target] = []

save_pickle('../results_mi/top_pathways_multinet_mi_diag.pickle', path_dict)