import pickle
import pandas as pd

def save_pickle(file, thing, message=False):
    pickle_out = open(file, "wb")
    pickle.dump(thing, pickle_out)
    pickle_out.close()
    if message:
        print(f'new pickle created \'{file}\'')

levels = ['GENETIC', 'MOLECULAR', 'PET', 'MRI', 'RISKFACTORS', 'PHENOTYPE']
node_dictionary = {}
for level in levels:
    df = pd.read_csv(f'../datasets/{level}.csv')
    df.set_index('RID', inplace=True)
    node_dictionary[level] = df.columns.tolist()
save_pickle('node_dictionary.pickle', node_dictionary)