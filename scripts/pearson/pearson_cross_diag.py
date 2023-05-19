import numpy as np
import pandas as pd
import time
from scipy.stats import pearsonr
from itertools import combinations_with_replacement, combinations
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

levels = ['GENETIC', 'MOLECULAR', 'PET', 'MRI', 'RISKFACTORS', 'PHENOTYPE']
diags = ['nl', 'mci', 'dementia']
combos = list(combinations(levels, 2))
print(f"Computing the pearson coefficient for multilevel networks.")
since = time.time()
for level1, level2 in combos:
    for diag in diags:
        print(level1, level2, diag)
        df1 = pd.read_csv(f'../datasets/{level1}_{diag}.csv')
        df1.set_index('RID', inplace=True)
        df2 = pd.read_csv(f'../datasets/{level2}_{diag}.csv')
        df2.set_index('RID', inplace=True)
        nvars1, nvars2 = df1.shape[1], df2.shape[1]

        df_joined = join_df(df1, df2, level1, level2)
        print('# subjects:  ', df_joined.shape[0]) # rows
        print('# variables: ', df_joined.shape[1]) # columns

        # only pairs that go between networks, pairs within networks should already be calculated
        tups1 = range(nvars1) # de 0 a nvar1-1
        tups2 = range(nvars1, nvars1 + nvars2) # de nvars1 a nvars1+nvars2-1
        tuples = [(x,y) for x in tups1 for y in tups2] # lista de tuplas

        coef, pval = coef_matrix_sklearn(df_joined, tuples=tuples)
        coef_array = coef[np.triu_indices(coef.shape[0])]
        pval_array = coef[np.triu_indices(coef.shape[0])]

        file = f'../results_pearson/pearson_array_{level1}_{level2}_{diag}.npy'
        print(file)
        print("")
        np.save(file, coef_array)
        file = f'../results_pearson/pval_array_{level1}_{level2}_{diag}.npy'
        print(file)
        print("")
        np.save(file, pval_array)
time_elapsed = time.time() - since
print(f'Total time: {time_elapsed:.1f}s')