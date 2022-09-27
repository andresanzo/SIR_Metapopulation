"""
@name  : Run_SIMS.py
@author: ANDRES
@date:   MARZO 19, 2020
@update: JUNIO 19, 2020

PROGRAMA PARA LEER LOS ARCHIVOS HDF5 CONTENIDOS EN LAS CARPETAS RUN_<NUMBER> Y GRAFICAR EL TOTAL DE INFECTADOS
EN LA RED METAPOBLACIONAL

Entradas:

    Graph_Name: El nombre del grafo definido en el archivo Settings.py. String
    M:          El número de parches a controlar
    Plot:       Variable boolena para decidir si se muestra o no la grafica. Bool
    
Para ejecutar este programa desde una consola poner: python3 Read_and_Plot.py  <Graph_Name>  <M>  <True or False>
donde Graph_Name es el nombre del grafo, M el número de nodos a controlar y True/False para mostrar la grafica o no.

Para compilar este programa desde una consola teclear: python3 Read_and_Plot.py <Nombre del grafo> <Numero de parches a controlar M> <Plot = True/False>

Para compilar este programa es necesario compilar primero el programa Rund_SIMS.py con las mismos parámetros 

"""

import h5py, sys, csv, pickle, glob, shutil, os
from subprocess    import call
from   tqdm                import tqdm
from   Settings            import *
from   SIR_Functions       import Save_Info
from   numpy               import linspace, arange, mean
import networkx            as nx
import matplotlib.pyplot   as plt
import pandas              as pd


# Formato del grafico
fontsize  = 15
linewidth = 2
labelsize = 15
alpha     = 0.6
formato   = 'png'

Graph_Name = sys.argv[1]
M          = int(sys.argv[2])
Plot       = sys.argv[3]

# Para graficar y guardar las series de tiempo de todas las corridas.
Plot_All   = False

# Arreglo para guardar la serie de los infectados cuando se controla el parche con el max I
Nun_Inf_ctr_1, Nun_Inf_ctr_2, Nun_Inf_ctr_5  = [], [], []

# Arreglo para guardar la serie de los recuperados cuando se controla el parche con el max I
Nun_R_ctr_1, Nun_R_ctr_2, Nun_R_ctr_5        = [], [], []

# Arreglo para guardar las propiedes de los nodos seleccionados con el criterio I_max
Nodes_prop_Imax = []

# Arreglo para guardar las propiedes de los nodos seleccionados aleatoriamente
Nodes_prop_Rand = []

# Arreglo para guardar el maximo de la suma de los infectados para cuando se controla aleatoriamente, con el índice y sin control
Infected_MAX    = []

Name_File = Graph_Name + '_M_' + str(M)

# Ciclo sobre cada run para guardar la informaci\'on
for run_i in tqdm(range(1,N_runs_ensamble_max)):  

        # Abriendo el archivo pkl
        infile         = open(Name_File + '/RUN_' + str(run_i) + '/'+ 'M_SELECTED_PATCH_RUN_' + str(run_i) + '.pkl','rb')
        new_dict       = pickle.load(infile)
        
        ctr_patch      = new_dict[0:M]
        random_patches = new_dict[M:]

        # ========================================================================================================================
        #    GRAFICA DE LOS INFECTADOS CONTROLANDO M PARCHEs CON EL M\'AXIMO PORCENTAJE DE INFECTADOS
        # ========================================================================================================================

        # Abrir el archivo Max_I_Run_<run_i>.h5py con la simulaci\'on controlando M  parches con el m\'aximo porcentaje de I
        file_1    = h5py.File(Name_File +'/RUN_' + str(run_i) + '/' + 'Max_I_Run_' + str(run_i) + '.hdf5' ,'r')["default"]
        data_1    = file_1[:int(px/hpx)]

        # Series de tiempo para los infectados (Controlando el parche con el m\'aximo porcentaje)
        I_set_1   = [data_1[:,3*vi+1] for vi in range(N)]

        # Series de tiempo para los recuperados (Controlando el parche con el m\'aximo porcentaje)
        R_set_1   = [data_1[:,3*vi+2] for vi in range(N)]
        
        # Realiza la suma de todas las series de tiempo de los infectados I_set_1 y lo gurada en el arreglo Nun_Inf_ctr_1
        I_zip_1   = zip(*[I_set_1[patchi] for patchi in range(N)])
        Nun_Inf_ctr_1.append( [sum(ei) for ei in I_zip_1] )

        # Realiza la suma de todas las series de tiempo de los recuperados R_set_1 y lo gurada en el arreglo Nun_R_ctr_1
        R_zip_1   = zip(*[R_set_1[patchi] for patchi in range(N)])
        Nun_R_ctr_1.append( [sum(ei) for ei in R_zip_1] )
        
        # ========================================================================================================================
        #    GRAFICA DE LOS INFECTADOS CONTROLANDO M PARCHES SELECCIONADOS AL AZAR
        # ========================================================================================================================

        # Abrir el archivo Random_Patches_Run_<run_i>.h5py
        file_2   = h5py.File(Name_File +'/RUN_' + str(run_i) + '/' + 'Random_Patches_Run_' + str(run_i) + '.hdf5' ,'r')["default"]
        data_2   = file_2[:int(px/hpx)]
    
        # Series de tiempo para los infectados (Controlando el parche con el m\'aximo porcentaje)
        I_set_2  = [data_2[:,3*vi+1] for vi in range(N)]

        # Series de tiempo para los recuperados (Controlando el parche con el m\'aximo porcentaje)
        R_set_2  = [data_2[:,3*vi+2] for vi in range(N)]
    
        # Realiza la suma de todas las series de tiempo de los infectados I_set_2 y lo gurada en el arreglo Nun_Inf_ctr_2
        I_zip_2  = zip(*[I_set_2[patchi] for patchi in range(N)])
        Nun_Inf_ctr_2.append( [sum(ei) for ei in I_zip_2] )

        # Realiza la suma de todas las series de tiempo de los recuperados R_set_2 y lo gurada en el arreglo Nun_R_ctr_2
        R_zip_2  = zip(*[R_set_2[patchi] for patchi in range(N)])
        Nun_R_ctr_2.append( [sum(ei) for ei in R_zip_2] )

        # ========================================================================================================================
        #    GRAFICA DE LOS INFECTADOS SIN CONTROL
        # ========================================================================================================================
    
        # Abrir el archivo Max_K_in_Run_<run_i>.h5py con la simulaci\'on controlando el parche con el m\'aximo porcentaje de I
        file_5   = h5py.File(Name_File +'/RUN_' + str(run_i) + '/' + 'Without_Control_Run_' + str(run_i) + '.hdf5' ,'r')["default"]
        data_5   = file_5[:int(px/hpx)]
    
        # Series de tiempo para los infectados (Controlando el parche con el m\'aximo porcentaje)
        I_set_5  = [data_5[:,3*vi+1] for vi in range(N)]

        # Series de tiempo para los infectados (Controlando el parche con el m\'aximo porcentaje)
        R_set_5  = [data_5[:,3*vi+2] for vi in range(N)]
    
        # Realiza la suma de todas las series de tiempo de los infectados I_set_2 y lo gurada en el arreglo Nun_Inf_ctr_5
        I_zip_5  = zip(*[I_set_5[patchi] for patchi in range(N)])
        Nun_Inf_ctr_5.append( [sum(ei)   for ei in I_zip_5] )

        # Realiza la suma de todas las series de tiempo de los recuperados R_set_2 y lo gurada en el arreglo Nun_R_ctr_5
        R_zip_5  = zip(*[R_set_5[patchi] for patchi in range(N)])
        Nun_R_ctr_5.append( [sum(ei)   for ei in R_zip_5] )
    
        # ========================================================================================================================
        #    Guardando la información en un archivo csv
        # ========================================================================================================================
        
        # Guardando la información sobre le máximo de la suma de los infectados
        Infected_MAX.append( [run_i,max(Nun_Inf_ctr_1[run_i-1]), max(Nun_Inf_ctr_2[run_i-1]), max(Nun_Inf_ctr_5[run_i-1]),
        					  max(Nun_R_ctr_1[run_i-1]), max(Nun_R_ctr_2[run_i-1]), max(Nun_R_ctr_5[run_i-1])	
        				     ] )
        
        # Abriendo el archivo pkl
        infile2   = open(Name_File +'/RUN_' + str(run_i) + '/' + 'INITIAL_CONDITIONS_RUN_' + str(run_i) +'.pkl','rb')
        new_dict2 = pickle.load(infile2)
        
        # La matriz de movilidad
        P                  = new_dict2[1]
        
        # La distribucion del total de la poblacion en cada parche
        Ni                  = new_dict2[4]
        
        # La matriz de movilidad
        G                  = new_dict2[5]
        
        # Diccionarios con las caracteristicas estructurales de los nodos
        
        # Compute the betweenness centrality for nodes in G: the fraction of number of shortests paths that pass through each node
        betweenness       = nx.betweenness_centrality(G)
        
        #Degree centrality for nodes (fraction of nodes connected to). Returns a dictionary of degree centrality values keyed by node.
        degree_centrality = nx.degree_centrality(G)
        
        # Closeness centrality for nodes (1/average distance to all nodes).
        closeness         = nx.closeness_centrality(G)
        
        # Guarda la información en arreglos de M vectores por cada run con las propiedades de los nodos. Este info será guardada en "TABLA_GRAPH" + graph_name +  ".csv"
        Nodes_prop_Imax.append( [  [vi,max(data_1[:,3*vi+1]),hpx*sum(data_1[:,3*vi+1]),Ni[vi],G.degree(vi),sum(P[:,vi]) - P[vi,vi],
                                    sum(P[vi,:]) -P[vi,vi],betweenness[vi],degree_centrality[vi],closeness[vi] ]   for vi in ctr_patch]    )
                                    
        Nodes_prop_Rand.append( [  [vi,max(data_2[:,3*vi+1]),hpx*sum(data_2[:,3*vi+1]),Ni[vi],G.degree(vi),sum(P[:,vi]) - P[vi,vi],
                                    sum(P[vi,:]) -P[vi,vi],betweenness[vi],degree_centrality[vi],closeness[vi] ]   for vi in random_patches]    )


# =================================================================================================================
#                  Rutina para escribir el máximo de los infectados global
# =================================================================================================================

with open( "MAX_Infected_" + graph_name +  ".csv",'w+',) as csvDataFile:

    filewriter = csv.writer(csvDataFile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    filewriter.writerow(['Run','Max_I_Index','Max_I_Random','Max_I_NoCtr','Effective','Max_R_Index','Max_R_Random','Max_R_NoCtr','Effective_R'])   # Headers
    
    for vi in Infected_MAX:
    
        if   vi[3] - vi[1] < vi[3] - vi[2]: effective_I = 1
        elif vi[3] - vi[1] > vi[3] - vi[2]: effective_I = 0

        if   vi[6] - vi[4] < vi[6] - vi[5]: effective_R = 1
        elif vi[6] - vi[4] > vi[6] - vi[5]: effective_R = 0
        
        filewriter.writerow([vi[0],vi[1],vi[2],vi[3],effective_I, vi[4],vi[5],vi[6], effective_R])
        
# =================================================================================================================
#           Rutina para escribir los encabezados del archivo Tablas  'Max_Infected','Total_Infected'
# =================================================================================================================

# Encabezados para el archivo "TABLA_GRAPH" + graph_name +  ".csv"
# Las Primeras columnas y filas para guardar la información sobre el nodo seleccionado con el algoritmo iterativo
first_columns = []
for vi in range(M):
    first_columns.extend(['Run','Selected Imax Patch_' + str(vi+1),'Max_Infected','Total_Infected','Ni', 'Node_degree','k_in','k_out','betweenness centrality','Degree centrality','Closeness centrality'   ])

# Las siguientes columnas para guardar la información sobre cada uno de los M parches seleccionados aleatoriamente
next_columns = []
for vi in range(M):
    next_columns.extend(['Run','Selected Random Patch_' + str(vi+1),'Max_Infected','Total_Infected','Ni', 'Node_degree','k_in','k_out','betweenness centrality','Degree centrality','Closeness centrality'   ])
    
# TABLA PARA GUARDAR LA INFORMACION SOBRE LAS CARACTERISTICAS ESTRUCTURALES DEL NODO A CONTROLAR
with open( "TABLA_GRAPH_" + graph_name +  ".csv",'w+',) as csvDataFile:

    filewriter = csv.writer(csvDataFile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)

    # Guardando la informaci\'on sobre las propiedades estructurales de los nodos seleccionados con el Max I_max (\'indice propuesto)
    filewriter.writerow(first_columns)   # Headers
    
    # # Guardando la informaci\'on sobre las propiedades estructurales de los M nodos seleccionados con el Max I_max
    for run,data_Imax in enumerate(Nodes_prop_Imax):
        
        column = []
        for i in range(M): column.extend([run+1]+data_Imax[i])
            
        filewriter.writerow(column)
    
    # Guardando la informaci\'on sobre las propiedades estructurales de los nodos seleccionados aleatoriamente
    filewriter.writerow(next_columns)  # Headers
    for run,data_random in enumerate(Nodes_prop_Rand):
    
        column = []
        for i in range(M): column.extend([run+1]+data_random[i])
    
        filewriter.writerow(column)

# =================================================================================================================
#         Fin de la rutina para escribir los encabezados del archivo Tablas  'Max_Infected','Total_Infected'
# =================================================================================================================


# Hace el promedio sobre cada realizaci\'on de los infectados controlando los M parches con el max I
Ctr1        = zip(*[Nun_Inf_ctr_1[runi] for runi in range(N_runs_ensamble_max-1)])
I_mean_ctr1 = [sum(ei)/(N_runs_ensamble_max-1)  for ei   in Ctr1 ]

# Hace el promedio sobre cada realizaci\'on de los infectados controlando aleatoriamente M parches
Ctr2        = zip(*[Nun_Inf_ctr_2[runi] for runi in range(N_runs_ensamble_max-1)])
I_mean_ctr2 = [sum(ei)/(N_runs_ensamble_max-1) for ei   in Ctr2 ]

# Hace el promedio sobre cada realizaci\'on de los infectados sin control
Ctr5        = zip(*[Nun_Inf_ctr_5[runi] for runi in range(N_runs_ensamble_max-1)])
I_mean_ctr5 = [sum(ei)/(N_runs_ensamble_max-1)  for ei   in Ctr5 ]

# ===============================================================================================================================================
#                                       SAVE INFO
# ===============================================================================================================================================
    
# Guardando los archivos csv en la carpeta Graph_Name
Save_Info().SaveFiles(Name_File,make_dir=False,file_format='*.csv')


# ===============================================================================================================================================
#                                   PLOT
# ===============================================================================================================================================

if Plot == 'True':

    plt.figure(10,figsize=(15,7))

    if M == 1:
        legend_for_Ctr_M_Patches =  str(M) + ' patch'
    elif M > 1:
        legend_for_Ctr_M_Patches =  str(M) + ' patches'

    plt.title(dict[graph_name])
    plt.plot(arange(0.0,px,hpx),I_mean_ctr1,color = 'red',linewidth=linewidth,alpha = alpha,label='Controlling ' + legend_for_Ctr_M_Patches + ' with criterium')
    plt.plot(arange(0.0,px,hpx),I_mean_ctr2,color = 'green',linewidth=linewidth,alpha = alpha,label='Controlling randomly ' + legend_for_Ctr_M_Patches)
    plt.plot(arange(0.0,px,hpx),I_mean_ctr5,color = 'black',linewidth=linewidth,alpha = alpha,label='Without control')

    plt.ylabel("$I_{i}(t)$",fontsize=fontsize)
    plt.xlabel("$time$",fontsize=fontsize)
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.legend(loc=1, bbox_to_anchor=(0.9,0.8),fontsize=10,ncol=1, fancybox=True, shadow=True)

    plt.show()


if Plot_All  == True:
    
    for runi in range(N_runs_ensamble_min,N_runs_ensamble_max):
    
        call(['python3','SIR_Paper_Figures.py',str(runi),Graph_Name])
        
    for ifile in glob.iglob(os.path.join('*.png')):
        shutil.move(ifile, Graph_Name + '_M_' + str(M))
    

