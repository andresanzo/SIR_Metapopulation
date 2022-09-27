# Generate_Initial_Conditions.py
# Fecha: Febrero 22, 2020
# Programa para generar el conjunto de condiciones iniciales para correr cada simulacion. 
# La salida del programa es un archivo pkl con los siguientes datos
#    -- El conjunto de condiciones iniciales de cada nodo (Nhi,Nvi y el parche infectado)
#    -- La matriz de movilidad P = (Pij)

import pickle, sys
from   numpy          import random
from   SIR_Functions  import Models
from   Settings       import *

# Numero de run
Run_Number = sys.argv[1]

# Condiciones iniciales para cada parche y el arreglo Ni con los valores del tamaño de la poblaci\'on en cada parche
Init_Conds,Ni      = Models().Init_Conditions(Nh_range=[Nh_min,Nh_max])

# La matriz de movilidad
P,G              = Models().P(N,graph_name,pij_random)

# Arreglo con la distribuci\'on de valores beta (tasa de contacto )de cada parche
betas = [random.uniform(beta_min,beta_max) for i in range(N)]

# Arreglo con la distribuci\'on de valores gamma (tasa de recuperaci\'on) de cada parche
gammas = [random.uniform(gamma_min,gamma_max) for i in range(N)]


# ========================================================
#	   Guardar la informacion con pickle
# ========================================================

Simulation_Data = [Init_Conds,P,betas,gammas,Ni,G]

outfile = open('INITIAL_CONDITIONS_RUN_' + Run_Number + '.pkl','wb')
pickle.dump(Simulation_Data,outfile)
outfile.close()

