# -*- coding: utf-8 -*-
"""
@name  : SIR_Controlled_M_Patches.py
@author: ANDRES
@date:   MARZO 18, 2020

 PROGRAMA PARA REALIZAR SIMULACIONES SELECCIONANDO ALEATORIAMENTE M NODOS EN LOS CUALES SE REALIZA EL CONTROL
 LA MATRIZ DE MOVILIDAD, LAS CONDICIONES INICIALES DEL SISTEMA Y LA DISTRIBUCI\'ON DE BETAS SON LEIDOS DEL ARCHIVO INITIAL_CONDITIONS_RUN_<NUM>.pkl
 ENTRADAS DEL PROGRAMA
    -- EL CONJUNTO DE NODOS A CONTROLAR
    -- EL NUMERO DE CORRIDA
    -- EL NOMBRE DEL ARCHIVO DE SALIDA
    
 LA SALIDA DEL SISTEMA
    -- UN ARCHIVO EN FORMATO HDF5 CON LAS SERIES DE TIEMPO DE LOS N PARCHES
"""

import sys,pickle,h5py
from   Settings            import *
from   SIR_Functions       import Models
from   scipy.integrate     import odeint
from   numpy               import random,array

# ====================================================================
#          ENTRADAS DEL PROGRAMA Y CONDICIONES INICIALES
# ====================================================================

# Arreglo con la lista de nodos a controlar
Controlled_Patches  = Models().StringArray_to_FloatArray(sys.argv[1])

# El n\'umero de corrida de la simulaci\'on
Run_Number          = sys.argv[2]

# EL nombre del archivo de salida
fname               = sys.argv[3]

# Si se desea un protocolo de control o no. Las opciones son: 'Sin_control', 'Control_Beta', 'Control_Pij_in', 'Control_Pij_out'
control_protocol    = sys.argv[4]

# Abriendo el archivo pkl
infile   = open('INITIAL_CONDITIONS_RUN_' + Run_Number +'.pkl','rb')
new_dict = pickle.load(infile)

# Condiciones iniciales
Initial_Conditions   = new_dict[0]

# La matriz de movilidad
P                    = new_dict[1]

# Los valores beta para cada parche
betas                = new_dict[2]

# Los valores gamma para cada parche
gammas               = new_dict[3]

# Los valores de la poblaci\'on total en cada parche
Ni                   =  new_dict[4]

# ================================================================================================================
#  RESOLVIENDO EL SISTEMA DE ECUACIONES CON EL CORRESPONDIENTE PROTOCOLO DE CONTROL Y GUARDANDO LA INFORMACI\'ON
# ================================================================================================================

if  control_protocol == 'Sin_control':

    X   = odeint( Models().SIR, Initial_Conditions,time,args=(P,betas,gammas,Ni) )
    
    #X_last = X[-1]
    #print('='*100)
    #print( [X_last[3*vi+2] for vi in range(N)])
    #print('='*100,'\n')
    
if  control_protocol == 'Control_Beta':
    
    for vi in Controlled_Patches: betas[vi] = (1-uc)*betas[vi]

    X   = odeint( Models().SIR, Initial_Conditions,time,args=(P,betas,gammas,Ni) )

if  control_protocol == 'Control_Pij_out':

    # Se realiza una copia de la matriz de movilidad
    P_new = P.copy()
    
    # El elemento de la diagonal de los  parches a controlar se hace uno y el resto de los elementos en las otras columnas se hacen cero. Se hace que el arreglo Controlled_Patches solo tiene un elemento.
    
    for vi in Controlled_Patches: # La etiqueta del parche a controlar
    
        # Todos las columnas de vi se reducen a un rango de valores pequeño
        P_new[vi,0:N]  = random.uniform(0.05,0.09)
    
        # El elemento de la diagonal del nodo a controlar corresponde a la suma del resto de las columnas
        P_new[vi,vi]   = sum( [P_new[vi,vj] for vj in range(N) if vj != vi ] )
    
    X   = odeint( Models().SIR, Initial_Conditions,time,args=(P_new,betas,gammas,Ni) )
    
if  control_protocol == 'Control_Pij_in':

    # Se realiza una copia de la matriz de movilidad
    P_new = P.copy()

    # El elemento de la diagonal del parche a controlar se hace uno y el resto de los elementos en las otras columnas se hacen cero. Se hace que el arreglo Controlled_Patches solo tiene un elemento.
    for vi in Controlled_Patches: # La etiqueta del parche a controlar
        for k in range(N):
    
            if k != vi:
        
                if P[k,vi] != 0:
        
                    # La movilidad hacia adentro del parche vi es reducida
                    P_new[k,vi] = random.uniform(0.05,0.07)
        
                    # Se modifia el elemento de la diagonal correspondiente al parche k
                    # El valor previo de su entrada en la diagonal mas lo que se le quito de movilidad hacia el parche seleccionado
                    P_new[k,k] = P[k,k] + (P[k,vi] - P_new[k,vi])
        
    X   = odeint( Models().SIR, Initial_Conditions,time,args=(P_new,betas,gammas,Ni) )
    
# ====================================================================
#               GUARDANDO LA INFORMACI\´ON
# ====================================================================

with h5py.File(fname+'.hdf5', 'w') as fh5:
    fh5.create_dataset("default",data=X)
    fh5.close()
