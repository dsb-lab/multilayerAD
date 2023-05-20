# Paper

Este repositorio contiene el código relacionado con el análisis y la simulación de redes para nuestro artículo de investigación.

## Creación del diccionario de nodos
- `node_dict.py`: Crea el archivo `node_dictionary.pickle`, donde las claves son las capas y los valores son los nodos.

## Mutual Information
### Redes Individuales
- `mi_ind.py`: Calcula la MI de las redes individuales y realiza el test de permutación.
- `mi_ind_net.py`: Crea las redes individuales a partir de los archivos generados por `mi_ind.py`.
### Redes Cross
- `mi_cross.py`: Calcula la MI de las redes cross y realiza el test de permutación.
- `mi_cross_net.py`: Crea las redes cross a partir de los archivos generados por `mi_cross.py`.
### Redes Individuales para nl, mci y dementia
- `mi_ind_diag.py`: Calcula la MI de las redes individuales para nl, mci y dementia y realiza el test de permutación.
- `mi_ind_diag_net.py`: Crea las redes individuales para nl, mci y dementia a partir de los archivos generados por `mi_ind_diag.py`.
### Redes Cross para nl, mci y dementia
- `mi_cross_diag.py`: Calcula la MI de las redes cross para nl, mci y dementia y realiza el test de permutación.
- `mi_cross_diag_net.py`: Crea las redes cross para nl, mci y dementia a partir de los archivos generados por `mi_cross_diag.py`.
### Análisis de la Red Global y Multicapa
- `bignet_analysis.ipynb`: Realiza el análisis de las redes individuales y crea la red global y multicapa.
### Cálculo de los 10 Shortest Paths
- `make_paths.py`: Calcula los 10 shortest paths para cada nodo inicial y final en la red global o multicapa.
- `make_paths_diag.py`: Calcula los 10 shortest paths para cada nodo inicial y final en las redes de nl, mci y dementia.
### Análisis de los Shortest Paths y Plots
- `paths_analysis.ipynb`: Realiza el análisis de los shortest paths y genera gráficas.
### Análisis de Conectividad y Matrices de Densidad
- `cross_analysis.ipynb`: Calcula las matrices de densidad (conectividad).


### PEARSON

`pearson_ind.py` -> calcula el pearson de las redes individuales

`pearson_ind_net.py` -> crea las redes individuales a partir de los files de `pearson_ind.py`

`pearson_cross.py` -> calcula el pearson de las redes cross

`pearson_cross_net.py` -> crea las redes cross a partir e los files de `pearson_cross.py`

`pearson_ind_diag.py` -> calcula el pearson de las redes individuales para nl, mci, dementia

`pearson_ind_net_diag.py` -> crea las redes individuales para nl, mci, dementia a partir de los files de `pearson_ind_diag.py`

`pearson_cross_diag.py` -> calcula el pearson de las redes cross para nl, mci, dementia

`pearson_cross_net_diag.py` -> crea las redes cross para nl, mci, dementia a partir e los files de `pearson_cross_diag.py`

`bignet_analysis_pearson.ipynb` -> analisis de las redes individuales y crea la red global y multicapa, también crea la red global para nl, mci, dementia


### SIMULACIONES BOOLEANAS

`boolean_simulations.py` -> Hace las simulaciones Booleanas, se puede tunear el noise

`boolean_simulations_diag.py` -> Hace las simulaciones Booleanas para nl, mci y dementia, se puede tunear el noise

`cross_corr.py` -> Calcula el cross correlation coefficient (maximum value of the temporal cross correlation across all lag times)

`cross_corr_diag.py` -> Calcula el cross correlation coefficient para nl, mci, dementia

`noise_analysis.ipynb` -> Analiza el efecto del noise en el cross correlation coefficient (falta representar la std)

`paths_cross_corr.py` -> Crea la network a partir del cross correlation coeficient y calcula los top shortest pathways

`paths_cross_corr_diag.py` -> Crea la network a partir del cross correlation coeficient y calcula los top shortest pathways para nl, mci, dementia

`permut_bool.py` -> Hace el negative control, permutando los edges de las redes y calculando las simulaciones booleanas y los top paths en cada permutation.

`analysis_paths.ipynb` -> Análisis de los paths obtenidos en el análisis dinámico. En los negative controls: contador de las veces que aparecen los paths después de hacer las permutaciones
