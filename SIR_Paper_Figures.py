# -*- coding: utf-8 -*-
"""
@name  : SIR_Fig1.py
@author: ANDRES
@date:   MAYO 25, 2020
@update: JUNIO 2, 2020

PROGRAMA PARA REALIZAR SIMULACIONES Y OBSERVAR LA DINÁMICA LOCAL DE LOS PARCHES CUANDO SE REALIZA LOS DISTINTOS PROTOCOLOS DE CONTROL.
EN ESTE PROGRAMA NO SE REALIZA LA SELECCIÓN DE FORMA ALEATORIA DE LOS PARCHES, SOLO SIGUIENDO CUALQUIERA DE LOS ÍNDICES (ABSOLUTO O RELATIVO).

EL PROGRAMA REALIZA LOS SIGUIENTES PASOS POR CADA SIMULACI\'ON

STEP 0:
        -- GENERA LA RED DE MOVILIDAD CON LA ESTRUTURA DEFINIDA EN EL ARCHIVO Settings.py
        -- GENERA LAS CONDICIONES INICIALES
        -- GENERA LA DISTRIBUCI\'ON DE BETAS
        -- TODO LO ANTERIOR LO REALIZA EJECUTANDO EL ARCHIVO Generate_Initial_Conditions.py  <Numero de Run=0>
        -- GENERA LA CARPETA RUN_<Numero de Run> EN DONDE SE GUARDAR\'AN LOS RESULTADOS
        
STEP 1:
        -- RESOLVER LAS ECUACIONES EN DIFERENCIAS Y ENCONTRAR EL PARCHE CON EL MAXIMO PORCENTAJE DE INFECTADOS
        
STEP 2:
        -- REALIZAR LA SIMULACI\'ON CONTROLANDO EL PAR\'AMETRO BETA DEL PARCHE CON EL MAXIMO PORCENTAJE DE INFECTADOS
        
STEP 3:
        -- REALIZAR LA SIMULACI\'ON CONTROLANDO EL PAR\'AMETRO Pji (EL FLUJO DE ENTRADA) DEL PARCHE CON EL MAXIMO PORCENTAJE DE INFECTADOS
        
STEP 4:
        -- REALIZAR LA SIMULACI\'ON CONTROLANDO EL PAR\'AMETRO Pij (EL FLUJO DE SALIDA) DEL PARCHE CON EL MAXIMO PORCENTAJE DE INFECTADOS
               
STEP 5:
     -- REALIZAR SIMULACI\'ON SIN CONTROL

STEP 7:
      -- GUARDAR LA INFORMACI\'ON EN LA CARPETA RUN_<run_i>
      
STEP 8:
      -- HACE LA GRÁFICA DE LOS INFECTADOS DE LOS PASOS 2 AL 5.
      
    
 ENTRADAS DEL SISTEMA
    -- NOMBRE DEL GRAFO

 LA SALIDA DEL SISTEMA
    -- UNA CARPERTA POR CADA CORRIDA CON LAS SERIES DE TIEMPO DE LAS SIMULACIONES
    -- LAS GRAFICAS DE LOS INFECTADOS DE LOS PASOS 2 AL 5
    
Para compilar este programa desde una consala teclear  python3 SIR_Fig4.py

"""

import h5py, pickle, sys
from Settings             import *
from numpy                import linspace, argmax
import matplotlib.pyplot   as plt

def num_format(number):
    return ("{:,}".format(round(number,2)))

# =================================================================================================
#                               SETTINGS
# =================================================================================================

# El número de corrida que se va a graficar
Run_Number  = int(sys.argv[1])

# El nombre del grafo
Graph_type = sys.argv[2]

# El numero de parches seleccionados
M_select    = int(1)

# =================================================================================================
#                               OPEN H5 FILES
# =================================================================================================

# El nombre de la carpeta con las archivos a graficar
File_Name      = Graph_type + '_M_' + str(M_select) + '/RUN_' + str(Run_Number)

# Abriendo el archivo pkl con la información de cuales fueron los nodos controlados
infile         = open( File_Name + '/M_SELECTED_PATCH_RUN_' + str(Run_Number) + '.pkl','rb')
new_dict       = pickle.load(infile)

# El conjunto de parches seleccionados con el índice propuesto
ctr_patch      = new_dict[0:M_select]
selected_patch = ctr_patch[0]

# El conjunto de parches seleccionados aleatoriamente  --> En esta versión considero que solo un parche es seleccionado
selected_rand  = new_dict[1]

# Abriendo el archivo hdf5 para la simulacion sin control
file_no_Ctr    = h5py.File(File_Name + '/Without_Control_Run_' +  str(Run_Number) +'.hdf5','r')["default"]
X_no_Ctr       = file_no_Ctr[:int(px/hpx)]

# Abriendo el archivo hdf5 para la simulacion controlando el parche con alguno de los índices (relativo o absoluto) - Protocolo 1
file_Ctr_Max_I = h5py.File(File_Name + '/Max_I_Run_' +  str(Run_Number) +'.hdf5','r')["default"]
X_Ctr_Max_I    = file_Ctr_Max_I[:int(px/hpx)]

# Abriendo el archivo hdf5 para la simulacion controlando aleatoriamente el parche.
file_Ctr_Random = h5py.File(File_Name + '/Random_Patches_Run_' +  str(Run_Number) +'.hdf5','r')["default"]
X_Ctr_Random    = file_Ctr_Random[:int(px/hpx)]

# =================================================================================================
#                               PLOT
# =================================================================================================

# Formato del grafico
fontsize  = 15
linewidth = 1.
labelsize = 12
alpha     = 0.6
formato   = 'png'

cmap      = plt.cm.binary                  #.gray #gist_rainbow #Spectral #hsv #jet
colors    = cmap(linspace(0,1.0,N))

plt.figure(30,figsize=(30,10)) #dpi=30
#ax        = fig.add_subplot(111, autoscale_on=False, xlim=(-1, 5), ylim=(-3, 5))

# =========================================================================
#                   SUB-PLOT 1
# =========================================================================

plt.subplot(3,3,1)
plt.title('a) Susceptibles',fontsize=15)

for vi in range(N):
    plt.plot(time, X_no_Ctr[:,3*vi], color=colors[vi],linewidth=2.,alpha=0.5)
    
# To plot the selected patch with color red
for select_i in ctr_patch:
    plt.plot(time, X_no_Ctr[:,3*select_i], color = 'red',  linewidth=2.5,alpha=1.0)
    
# To plot the selected patch selected randomly with color green
plt.plot(time, X_no_Ctr[:,3*selected_rand],color = 'green',  linewidth=2.5,alpha=1.0)
    
plt.ylabel(r'Total population',fontsize=fontsize)
plt.tick_params(axis='y',labelsize=labelsize)
plt.tick_params(axis='x',labelsize=labelsize)
plt.xlim([0,20])

# =========================================================================
#                   SUB-PLOT 2
# =========================================================================

plt.subplot(3,3,2)
plt.title('b) Infected',fontsize=15)

for vi in range(N):
    plt.plot(time, X_no_Ctr[:,3*vi+1], color=colors[vi],linewidth=2.,alpha=0.5)
    
# To plot the selected patch with color red
for select_i in ctr_patch:
    plt.plot(time, X_no_Ctr[:,3*select_i+1], color = 'red',linewidth=2.5,alpha=1.0)

# Tiempo para llegar al máximo de los infectados
time_max_select_i = argmax(X_no_Ctr[:,3*selected_patch+1])*hpx

plt.scatter(time_max_select_i,max(X_no_Ctr[:,3*selected_patch+1]),marker='o',color='red',s=30)
plt.scatter(time_max_select_i,max(X_no_Ctr[:,3*selected_rand+1]),marker='o',color='green',s=30)

plt.annotate( num_format( max(X_no_Ctr[:,3*selected_patch+1]) ), xy=(time_max_select_i, max(X_no_Ctr[:,3*selected_patch+1])),  xycoords='data',
xytext=(0.9, 0.95), textcoords='axes fraction',fontsize=12,
horizontalalignment='right', verticalalignment='top',
arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})
    
# To plot the selected patch selected randomly with color green
plt.plot(time, X_no_Ctr[:,3*selected_rand+1],color = 'green',  linewidth=2.5,alpha=1.0)

plt.annotate( num_format( max(X_no_Ctr[:,3*selected_rand+1]) ), xy=(time_max_select_i, max(X_no_Ctr[:,3*selected_rand+1])),  xycoords='data',
xytext=(0.9, 0.75), textcoords='axes fraction',fontsize=12,
horizontalalignment='right', verticalalignment='top',
arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})
    
plt.tick_params(axis='y',labelsize=labelsize)
plt.tick_params(axis='x',labelsize=labelsize)
plt.xlim([0,20])

# =========================================================================
#                   SUB-PLOT 3
# =========================================================================

plt.subplot(3,3,3)
plt.title('c) Recovery',fontsize=15)

for vi in range(N):
    plt.plot(time, X_no_Ctr[:,3*vi+2], color=colors[vi],linewidth=2.,alpha=0.5)
    
# To plot the selected patch with color red
for select_i in ctr_patch:
    plt.plot(time, X_no_Ctr[:,3*select_i+2], color = 'red',linewidth=2.5,alpha=1.0)
    

plt.annotate( num_format( max(X_no_Ctr[:,3*selected_patch+2]) ), xy=(17.5, max(X_no_Ctr[:,3*selected_patch+2])),  xycoords='data',
xytext=(0.5, 0.95), textcoords='axes fraction',fontsize=12,
horizontalalignment='right', verticalalignment='top',
arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

# To plot the selected patch selected randomly with color green
plt.plot(time, X_no_Ctr[:,3*selected_rand+2],color = 'green',  linewidth=2.5,alpha=1.0)

plt.annotate( num_format( max(X_no_Ctr[:,3*selected_rand+2]) ), xy=(17.5, max(X_no_Ctr[:,3*selected_rand+2])),  xycoords='data',
xytext=(0.5, 0.5), textcoords='axes fraction',fontsize=12,
horizontalalignment='right', verticalalignment='top',
arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

#plt.ylabel(r'Number of infected',fontsize=fontsize)
plt.tick_params(axis='y',labelsize=labelsize)
plt.tick_params(axis='x',labelsize=labelsize)

plt.xlim([0,20])

# =========================================================================
#                   SUB-PLOT 4
# =========================================================================

plt.subplot(3,3,4)
#plt.title('d) Control protocol 3',fontsize=15)

for vi in range(N):
    plt.plot(time, X_Ctr_Max_I[:,3*vi], color=colors[vi],linewidth=2.,alpha=0.5)
    
# To plot the selected patch with color red
for select_i in ctr_patch:
    plt.plot(time, X_Ctr_Max_I[:,3*select_i], color = 'red',linewidth=2.5,alpha=1.0)
    
# To plot the selected patch selected randomly with color green
plt.plot(time, X_Ctr_Max_I[:,3*selected_rand],color = 'green',  linewidth=2.5,alpha=1.0)
   
plt.ylabel(r'Total population',fontsize=fontsize)
plt.tick_params(axis='y',labelsize=labelsize)
plt.tick_params(axis='x',labelsize=labelsize)

plt.xlim([0,20])


# =========================================================================
#                   SUB-PLOT 5
# =========================================================================

plt.subplot(3,3,5)
#plt.title('d) Control protocol 3',fontsize=15)

for vi in range(N):
    plt.plot(time, X_Ctr_Max_I[:,3*vi+1], color=colors[vi],linewidth=2.,alpha=0.5)
    
# To plot the selected patch with color red
for select_i in ctr_patch:
    plt.plot(time, X_Ctr_Max_I[:,3*select_i+1], color = 'red',linewidth=2.5,alpha=1.0)

# Tiempo para llegar al máximo de los infectados
time_max_select_max = argmax(X_Ctr_Max_I[:,3*selected_patch+1])*hpx

plt.scatter(time_max_select_max,max(X_Ctr_Max_I[:,3*selected_patch+1]),marker='o',color='red',s=30)
plt.scatter(time_max_select_max,max(X_Ctr_Max_I[:,3*selected_rand+1]),marker='o',color='green',s=30)

plt.annotate( num_format( max(X_Ctr_Max_I[:,3*selected_patch+1]) ), xy=(time_max_select_max, max(X_Ctr_Max_I[:,3*selected_patch+1])),  xycoords='data',
xytext=(0.9, 0.95), textcoords='axes fraction',fontsize=12,
horizontalalignment='right', verticalalignment='top',
arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

# To plot the selected patch selected randomly with color green
plt.plot(time, X_Ctr_Max_I[:,3*selected_rand+1],color = 'green',  linewidth=2.5,alpha=1.0)

plt.annotate( num_format( max(X_Ctr_Max_I[:,3*selected_rand+1]) ), xy=(time_max_select_max, max(X_Ctr_Max_I[:,3*selected_rand+1])),  xycoords='data',
xytext=(0.9, 0.75), textcoords='axes fraction',fontsize=12,
horizontalalignment='right', verticalalignment='top',
arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

plt.tick_params(axis='y',labelsize=labelsize)
plt.tick_params(axis='x',labelsize=labelsize)

plt.xlim([0,20])

# =========================================================================
#                   SUB-PLOT 6
# =========================================================================

plt.subplot(3,3,6)
#plt.title('d) Control protocol 3',fontsize=15)

for vi in range(N):
    plt.plot(time, X_Ctr_Max_I[:,3*vi+2], color=colors[vi],linewidth=2.,alpha=0.5)
    
# To plot the selected patch with color red
for select_i in ctr_patch:
    plt.plot(time, X_Ctr_Max_I[:,3*select_i+2], color = 'red',linewidth=2.5,alpha=1.0)
    
plt.annotate( num_format( max(X_Ctr_Max_I[:,3*selected_patch+2]) ), xy=(17.5, max(X_Ctr_Max_I[:,3*selected_patch+2])),  xycoords='data',
xytext=(0.5, 0.95), textcoords='axes fraction',fontsize=12,
horizontalalignment='right', verticalalignment='top',
arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})
    
# To plot the selected patch selected randomly with color green
plt.plot(time, X_Ctr_Max_I[:,3*selected_rand+2],color = 'green',  linewidth=2.5,alpha=1.0)

plt.annotate( num_format( max(X_Ctr_Max_I[:,3*selected_rand+2]) ), xy=(17.5, max(X_Ctr_Max_I[:,3*selected_rand+2])),  xycoords='data',
xytext=(0.5, 0.5), textcoords='axes fraction',fontsize=12,
horizontalalignment='right', verticalalignment='top',
arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})
    
plt.tick_params(axis='y',labelsize=labelsize)
plt.tick_params(axis='x',labelsize=labelsize)

plt.xlim([0,20])

# =========================================================================
#                   SUB-PLOT 7
# =========================================================================

plt.subplot(3,3,7)

for vi in range(N):
    plt.plot(time, X_Ctr_Random[:,3*vi], color=colors[vi],linewidth=2.,alpha=0.5)
    
# To plot the selected patch with color red
for select_i in ctr_patch:
    plt.plot(time, X_Ctr_Random[:,3*select_i], color = 'red',  linewidth=2.5,alpha=1.0)
    
# To plot the selected patch selected randomly with color green
plt.plot(time, X_Ctr_Random[:,3*selected_rand],color = 'green',  linewidth=2.5,alpha=1.0)
    
plt.ylabel(r'Total population',fontsize=fontsize)
plt.xlabel(r'time (days)',fontsize=fontsize)
plt.tick_params(axis='y',labelsize=labelsize)
plt.tick_params(axis='x',labelsize=labelsize)
plt.xlim([0,20])

# =========================================================================
#                   SUB-PLOT 8
# =========================================================================

plt.subplot(3,3,8)

for vi in range(N):
    plt.plot(time, X_Ctr_Random[:,3*vi+1], color=colors[vi],linewidth=2.,alpha=0.5)
    
# To plot the selected patch with color red
for select_i in ctr_patch:
    plt.plot(time, X_Ctr_Random[:,3*select_i+1], color = 'red',linewidth=2.5,alpha=1.0)

# Tiempo para llegar al máximo de los infectados
time_max_select_rand = argmax(X_Ctr_Random[:,3*selected_patch+1])*hpx

plt.scatter(time_max_select_rand,max(X_Ctr_Random[:,3*selected_patch+1]),marker='o',color='red',s=30)
plt.scatter(time_max_select_rand,max(X_Ctr_Random[:,3*selected_rand+1]),marker='o',color='green',s=30)

plt.annotate( num_format( max(X_Ctr_Random[:,3*selected_patch+1]) ), xy=(time_max_select_rand, max(X_Ctr_Random[:,3*selected_patch+1])),  xycoords='data',
xytext=(0.9, 0.95), textcoords='axes fraction',fontsize=12,
horizontalalignment='right', verticalalignment='top',
arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})
    
# To plot the selected patch selected randomly with color green
plt.plot(time, X_Ctr_Random[:,3*selected_rand+1],color = 'green',  linewidth=2.5,alpha=1.0)

plt.annotate( num_format( max(X_Ctr_Random[:,3*selected_rand+1]) ), xy=(time_max_select_rand, max(X_Ctr_Random[:,3*selected_rand+1])),  xycoords='data',
xytext=(0.9, 0.75), textcoords='axes fraction',fontsize=12,
horizontalalignment='right', verticalalignment='top',
arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})
    
plt.xlabel(r'time (days)',fontsize=fontsize)
plt.tick_params(axis='y',labelsize=labelsize)
plt.tick_params(axis='x',labelsize=labelsize)
plt.xlim([0,20])

# =========================================================================
#                   SUB-PLOT 9
# =========================================================================

plt.subplot(3,3,9)

for vi in range(N):
    plt.plot(time, X_Ctr_Random[:,3*vi+2], color=colors[vi],linewidth=2.,alpha=0.5)
    
# To plot the selected patch with color red
for select_i in ctr_patch:
    plt.plot(time, X_Ctr_Random[:,3*select_i+2], color = 'red',linewidth=2.5,alpha=1.0)
    
plt.annotate( num_format( max(X_Ctr_Random[:,3*selected_patch+2]) ), xy=(17.5, max(X_Ctr_Random[:,3*selected_patch+2])),  xycoords='data',
xytext=(0.5, 0.95), textcoords='axes fraction',fontsize=12,
horizontalalignment='right', verticalalignment='top',
arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

# To plot the selected patch selected randomly with color green
plt.plot(time, X_Ctr_Random[:,3*selected_rand+2],color = 'green',  linewidth=2.5,alpha=1.0)

plt.annotate( num_format( max(X_Ctr_Random[:,3*selected_rand+2]) ), xy=(17.5, max(X_Ctr_Random[:,3*selected_rand+2])),  xycoords='data',
xytext=(0.5, 0.5), textcoords='axes fraction',fontsize=12,
horizontalalignment='right', verticalalignment='top',
arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

#plt.ylabel(r'Number of infected',fontsize=fontsize)
plt.xlabel(r'time (days)',fontsize=fontsize)
plt.tick_params(axis='y',labelsize=labelsize)
plt.tick_params(axis='x',labelsize=labelsize)

plt.xlim([0,20])

plt.savefig('SIR_RUN_'+ str(Run_Number) + '.png')

del time_max_select_i,time_max_select_max, time_max_select_rand

#plt.show()

# =================================================================================================
#                              PRINT INFO
# =================================================================================================

#from tabulate import tabulate

#print('---------------->',max(X_no_Ctr[:,3*selected_patch+1]))
#print('---------------->',max(X_Ctr_Max_I[:,3*selected_patch+1]))

##print(tabulate([['1', abs(max(X_no_Ctr[:,3*selected_patch+1]) - max(X_Ctr_Max_I[:,3*selected_patch+1]) ) ],
##                ['2', abs(max(X_no_Ctr[:,3*selected_patch+1]) - max(X_Ctr_K_in[:,3*selected_patch+1])  ) ]
##                ['3', abs(max(X_no_Ctr[:,3*selected_patch+1]) - max(X_Ctr_K_out[:,3*selected_patch+1]) ) ]],
##                headers=['Protocol', 'Differences among maxima']))

#print('='*100)
#print( 'Protocolo_1: ', abs(max(X_no_Ctr[:,3*selected_patch+1]) - max(X_Ctr_Max_I[:,3*selected_patch+1]) ) )

##print( 'Protocolo_2:   ', abs(max(X_no_Ctr[:,3*selected_patch+1]) - max(X_Ctr_K_in[:,3*selected_patch+1])  ) )
##print( 'Protocolo_3:  ', abs(max(X_no_Ctr[:,3*selected_patch+1]) - max(X_Ctr_K_out[:,3*selected_patch+1]) ) )
##print('='*100)

#plt.show()


