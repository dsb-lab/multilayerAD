# paper

node_dict.py -> crea el arachivo node_dictionary.pickle donde las keys son cada capa y los values son los nodos

MUTUAL INFORMATION

mi_ind.py -> calcula la MI de las redes individuales y hace el permutation test
mi_ind_net.py -> crea las redes individuales a partir de los files de mi_ind.py
mi_cross.py -> calcula la MI de las redes cross y hace el permutation test
mi_cross_net.py -> crea las redes cross a partir de los files de mi_cross.py
mi_ind_diag.py -> calcula la MI de las redes individuales para nl, mci, dementia y hace el permutation test
mi_ind_diag_net.py -> crea las redes individuales para nl, mci, dementia a partir de los files de mi_ind_diag.py
mi_cross_diag.py -> calcula la MI de las redes cross para nl, mci, dementia y hace el permutation test
mi_cross_diag_net.py -> crea las redes cross para nl, mci, dementia a partir de los files de mi_cross_diag.py

bignet_analysis.ipynb -> analisis de las redes individuales y crea la red global y multicapa

make_paths.py -> calcula los 10 shortest paths para cada nodo inicial y final (puede ser para la res global o la multicapa)
make_paths_diag.py -> igual que el anterior pero para nl, mci, dementia
paths_analysis.ipynb -> analisis de los shortest paths y plots

cross_analysis.ipynb -> calcula las matrices de densidad (conectividad)


PEARSON

pearson_ind.py -> calcula el pearson de las redes individuales
pearson_ind_net.py -> crea las redes individuales a partir de los files de pearson_ind.py
pearson_cross.py -> calcula el pearson de las redes cross
pearson_cross_net.py -> crea las redes cross a partir e los files de pearson_cross.py
pearson_ind_diag.py -> calcula el pearson de las redes individuales para nl, mci, dementia
pearson_ind_net_diag.py -> crea las redes individuales para nl, mci, dementia a partir de los files de pearson_ind_diag.py
pearson_cross_diag.py -> calcula el pearson de las redes cross para nl, mci, dementia
pearson_cross_net_diag.py -> crea las redes cross para nl, mci, dementia a partir e los files de pearson_cross_diag.py

bignet_analysis_pearson.ipynb -> analisis de las redes individuales y crea la red global y multicapa, también crea la red global para nl, mci, dementia


SIMULACIONES BOOLEANAS

boolean_simulations.py -> adaptado de Elena_net_paths. Hace las simulaciones Booleanas, se puede tunear el noise
boolean_simulations_diag.py -> Hace las simulaciones Booleanas para nl, mci y dementia, se puede tunear el noise
cross_corr.py -> adaptado de Elena_crosscor (primera parte). Calcula el cross correlation coefficient (maximum value of the temporal cross correlation across all lag times)
cross_corr_diag.py -> Calcula el cross correlation coefficient para nl, mci, dementia

noise_analysis.ipynb -> Analiza el efecto del noise en el cross correlation coefficient (falta representar la std)

paths_cross_corr.py -> adaptado de Elena_crosscor (segunda parte). Crea la network a partir del cross correlation coeficient y calcula los top shortest pathways
paths_cross_corr_diag.py -> Crea la network a partir del cross correlation coeficient y calcula los top shortest pathways para nl, mci, dementia

permut_bool.py -> adaptado de edge_swap.ipynb. Hace el negative control, permutando los edges de las redes y calculando las simulaciones booleanas y los top paths en cada permutation.

analysis_paths.ipynb -> Análisis de los paths obtenidos en el análisis dinámico. En los negative controls: contador de las veces que aparecen los paths después de hacer las permutaciones
