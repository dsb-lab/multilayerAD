import numpy as np
import pandas as pd
import time
from sklearn.metrics import mutual_info_score
from itertools import combinations_with_replacement
import warnings
warnings.filterwarnings("ignore")

def mi_calc_sklearn(z, bins=10):
    """
    calculate the mutual information between a pair of variables
    z is a 2d numpy array, where rows are variables, columns are subjects
    """
    nan_inds = np.where(np.isnan(z))[1] # find subject with NaN values
    z = np.delete(z, nan_inds, axis=1) # remove subjects with NaN values
    c_xy = np.histogram2d(z[0], z[1], bins=bins)[0]
    try:
        MI = mutual_info_score(None, None, contingency=c_xy) #mutual information
    except:
        MI = 0
    return MI

def coef_matrix_sklearn(df, bins=10, tuples=None, progress=True):
    """
    calculate the coefficient (mutual information or pearson) between all pairs in a pandas dataframe
    df       --> dataframe where rows are the subjects, columns are the variables
    bins     --> number of bins to use in MI calculation
    tuples   --> list of tuples cooresponding to which pairs to calculate MI
               (if not defined, all pairs are used)
    """
    n_vars = df.shape[1] # df dimension (rows x columns)
    if tuples == None:
        tuples = list(combinations_with_replacement(range(n_vars), 2))
        # no se repiten combinaciones (si tenemos (0,5) luego no tenemos (5,0))
    coef = np.zeros((n_vars, n_vars)) # define empty matrix (MI or pearson scores will go here)
    # tiene esta dimension porque es la MI entre cada variable con el resto
    for i, tupl in enumerate(tuples):
        if progress:
            print(f'pairs calculated: {i+1}/{len(tuples)}', end='\r')
        z = np.array(df.iloc[:, list(tupl)], dtype=float).T
        try:
            coef[tupl] = mi_calc_sklearn(z, bins=bins) # use function to calculate MI for pair
        except:
            coef[tupl] = 0
    coef = coef.T + np.triu(coef, k=1) # fill in the lower triangle of matrix
    return coef

def coef_perm_matrix(df, reps=1000, bins=10, tuples=None):
    nvars = df.shape[1]
    if tuples == None:
        tuples = list(combinations_with_replacement(range(nvars), 2))
    print('')
    print('# of variables:', nvars)
    print('# of pairs:', len(tuples))
    coef = np.zeros((reps, nvars, nvars)) # array 3d
    for rep in range(reps):
        for tupl in tuples:
            # por ej tupl=(3,4), list(tupl) da una lista con los dos números, así que z es un array con todas las rows
            # (pacientes) y las columnas 3 y 4 (las dos variables de las que vamos a calcular MI)
            z = np.array(df.iloc[:, list(tupl)], dtype=float).T # traspuesta para usarla de argumento en mi_calc_sklearn
            z = z[:, ~np.any(np.isnan(z), axis=0)]
            z = np.vstack((np.random.permutation(z[0]), np.random.permutation(z[1])))
            coef[rep, tupl[0], tupl[1]] = mi_calc_sklearn(z, bins=bins) # array 3d, por eso MI(rep, var1, var2)
        print(f'current rep: {rep+1}/{reps}', end='\r')
    return coef

levels = ['GENETIC', 'MOLECULAR', 'PET', 'MRI', 'RISKFACTORS', 'PHENOTYPE']
diags = ['nl', 'mci', 'dementia']
print("Computing the MI coefficient.")
since = time.time()
for level in levels:
    for diag in diags:
        print(level, diag)
        df = pd.read_csv(f'{level}_{diag}.csv')
        df.set_index('RID', inplace=True)
        print('# subjects:  ', df.shape[0]) # rows
        print('# variables: ', df.shape[1]) # columns
        coef = coef_matrix_sklearn(df, bins=10) # creates a square matrix of coef scores with shape (# variables, # variables)
        coef_array = coef[np.triu_indices(coef.shape[0])] # MI[np.triu_indices(rows)]
        file = f'{level}_{diag}_mi_array.npy'
        print(file)
        np.save(file, coef_array)

        coef = coef_perm_matrix(df, reps=1000, bins=10)
        coef_array = np.asarray([arr[np.triu_indices(arr.shape[0])] for arr in coef])
        file = f'{level}_{diag}_mi_perms.npy'
        print(file)
        np.save(file, coef_array)
time_elapsed = time.time() - since
print(f'Total time: {time_elapsed:.1f}s')