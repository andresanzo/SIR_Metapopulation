# Settings.py
# FECHA: ENERO 16, 2020
# Valores de los par\'ametros que se usar\'an en la simulaci\'on

# /Users/Golem/Documents/WORKING_PAPERS_2020/An_Iterative_Algoritm_Final_Size/PYTHON_INDICE_DE_INVASION_PAPER_UVENCIO/N_PATCHES_MAYO_2020

from numpy import arange

# NUMBER OF PATCHES
N                        = 15

# RANGO DE SIMULACIONES A REALIZAR

N_runs_ensamble_min      = 1

N_runs_ensamble_max      = 101

# El rango de valores para la densidad poblacional de cada parche
Nh_min, Nh_max           =  10000, 1000000

# El rango de valores para el par\'ametro beta
beta_min, beta_max       = 1.7, 2.5

# El rango de valores para el par\'ametro gamma
gamma_min, gamma_max     = 0.7,0.7

# NUMERO DE PARCHES A SELECCIONAR DE FORMA ALEATORIA
#M   = 2

# PARAMETRO DE CONTROL
uc = 0.5

# TOTAL DE PASOS PARA RESOLVER LAS ECUACIONE SIR
px = 30.0

# PASO DE INTEGRACION
hpx = 0.01

# RANGO DE TIEMPO
time = arange(0,px,hpx)

# SELECCIONA LA TOPOLOGIA DEL GRAFO. LAS OPCIONES SON: 'Fully', 'Ring', 'ER', 'NWS', 'BA'
graph_name      = 'BA'

# SELECCIONA LA DISTRIBUCION DE LOS PARAMETROS DE MOVILIDAD Pij. LAS OPCIONES SON: 'uniform', 'power', 'my_P'
pij_random      = 'uniform'

# PROTOCOLO PARA SELECCIONAR EL NODO A CONTROLAR: LAS OPCIONES SON 'relative' y 'absolute'
select_protocol = 'relative'

# SELECCIONA LA DISTRIBUCION DEL TAMAÑO DE LA POBLACION EN CADA PARCHE. LAS OPCIONES SON: 'uniform', 'power', 'my_Nh'
#Nh_random    = 'uniform'

# The form of select a patch where the infected host will be introduced. Option 'by_Nh' or 'by_fitness'
#select_patch = 'by_Nh'

dict ={'ER':'Erdos-Renyi', 'NWS': 'Newman-Watts-Strogatz', 'BA': 'Barabasi-Albert','Fully':'Fully connected network'}
