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

## Pearson
- `pearson_ind.py`: Calcula el coeficiente de correlación de Pearson para las redes individuales.
- `pearson_ind_net.py`: Crea las redes individuales a partir de los archivos generados por `pearson_ind.py`.
- `pearson_cross.py`: Calcula el coeficiente de correlación de Pearson para las redes cruzadas.
- `pearson_cross_net.py`: Crea las redes cruzadas a partir de los archivos generados por `pearson_cross.py`.
- `pearson_ind_diag.py`: Calcula el coeficiente de correlación de Pearson para las redes individuales de nl, mci y dementia.
- `pearson_ind_net_diag.py`: Crea las redes individuales de nl, mci y dementia a partir de los archivos generados por `pearson_ind_diag.py`.
- `pearson_cross_diag.py`: Calcula el coeficiente de correlación de Pearson para las redes cruzadas de nl, mci y dementia.
- `pearson_cross_net_diag.py`: Crea las redes cruzadas de nl, mci y dementia a partir de los archivos generados por `pearson_cross_diag.py`.
- `bignet_analysis_pearson.ipynb`: Realiza el análisis de las redes individuales y crea la red global y multicapa, también crea la red global para nl, mci y dementia.

## Simulaciones Booleanas
- `boolean_simulations.py`: Realiza simulaciones booleanas con la opción de ajustar el nivel de ruido.
- `boolean_simulations_diag.py`: Realiza simulaciones booleanas para nl, mci y dementia con la opción de ajustar el nivel de ruido.
- `cross_corr.py`: Calcula el coeficiente de correlación cruzada (maximum value of the temporal cross correlation across all lag times).
- `cross_corr_diag.py`: Calcula el coeficiente de correlación cruzada para nl, mci, dementia.
- `noise_analysis.ipynb`: Analiza el efecto del ruido en el coeficiente de correlación cruzada (falta representar la desviación estándar).
- `paths_cross_corr.py`: Crea la red a partir del coeficiente de correlación cruzada y calcula los mejores caminos más cortos (shortest paths).
- `paths_cross_corr_diag.py`: Crea la red a partir del coeficiente de correlación cruzada y calcula los mejores caminos más cortos para nl, mci, dementia.
- `permut_bool.py`: Realiza el control negativo, permutando las conexiones de las redes y calculando las simulaciones booleanas y los mejores caminos en cada permutación.
- `analysis_paths.ipynb`: Análisis de los caminos obtenidos en el análisis dinámico. Para los controles negativos: contador de las veces que aparecen los caminos después de realizar las permutaciones.
- `create_network_paths.ipynb`: A partir del CSV con los paths que han pasado el negative control crea las redes que represento en Cytoscape.