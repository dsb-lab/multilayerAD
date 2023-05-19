import numpy as np
import pandas as pd
import time
from scipy.stats import pearsonr
from itertools import combinations_with_replacement
import warnings
warnings.filterwarnings("ignore")

def coef_matrix_sklearn(df, tuples=None, progress=True):
    """
    calculate the coefficient (mutual information or pearson) between all pairs in a pandas dataframe
    df       --> dataframe where rows are the subjects, columns are the variables
    bins     --> number of bins to use in pearson calculation
    tuples   --> list of tuples cooresponding to which pairs to calculate pearson
               (if not defined, all pairs are used)
    """
    n_vars = df.shape[1] # df dimension (rows x columns)
    if tuples == None:
        tuples = list(combinations_with_replacement(range(n_vars), 2))
        # no se repiten combinaciones (si tenemos (0,5) luego no tenemos (5,0))
    coef = np.zeros((n_vars, n_vars))
    pval = np.zeros((n_vars, n_vars))
    for i, tupl in enumerate(tuples):
        if progress:
            print(f'pairs calculated: {i+1}/{len(tuples)}', end='\r')
        z = np.array(df.iloc[:, list(tupl)], dtype=float).T
        z = np.delete(z, np.where(np.isnan(z))[1], axis=1)  # delete rows with NaN indices
        try:
            coef[tupl], pval[tupl] = pearsonr(z[0], z[1])
        except:
            coef[tupl], pval[tupl] = 0, 1
    coef = coef.T + np.triu(coef, k=1) # fill in the lower triangle of matrix
    return coef, pval

levels = ['GENETIC', 'MOLECULAR', 'PET', 'MRI', 'RISKFACTORS', 'PHENOTYPE']
print(f"Computing the pearson coefficient.")
since = time.time()
for level in levels:
    print(level)
    df = pd.read_csv(f'../datasets/{level}.csv')
    df.set_index('RID', inplace=True)
    print('# subjects:  ', df.shape[0]) # rows
    print('# variables: ', df.shape[1]) # columns
    coef, pval = coef_matrix_sklearn(df) # creates a square matrix of coef scores with shape (# variables, # variables)
    coef_array = coef[np.triu_indices(coef.shape[0])]
    pval_array = pval[np.triu_indices(pval.shape[0])]
    file = f'../results_pearson/{level}_pearson_array.npy'
    print(file)
    np.save(file, coef_array)
    file = f'../results_pearson/{level}_pval_array.npy'
    print(file)
    np.save(file, pval_array)
time_elapsed = time.time() - since
print(f'Total time: {time_elapsed:.1f}s')