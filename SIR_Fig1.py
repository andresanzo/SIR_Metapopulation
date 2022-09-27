# -*- coding: utf-8 -*-
"""
@name  : SIR_Fig1.py
@author: ANDRES
@date:   ABRIL 15, 2020

 PROGRAMA PARA REALIZAR SIMULACIONES SELECCIONANDO ALEATORIAMENTE M NODOS EN LOS CUALES SE REALIZA EL CONTROL
 LA MATRIZ DE MOVILIDAD, LAS CONDICIONES INICIALES DEL SISTEMA Y LA DISTRIBUCI\'ON DE BETAS SON LEIDOS DEL ARCHIVO INITIAL_CONDITIONS_RUN_<NUM>.pkl
 ENTRADAS DEL PROGRAMA
    -- EL CONJUNTO DE NODOS A CONTROLAR
    -- EL NUMERO DE CORRIDA
    -- EL NOMBRE DEL ARCHIVO DE SALIDA
    
 LA SALIDA DEL SISTEMA
    -- UN ARCHIVO EN FORMATO HDF5 CON LAS SERIES DE TIEMPO DE LOS N PARCHES

"""

from Settings            import *
from SIR_Functions       import Models
from scipy.integrate     import odeint
from numpy               import random,array,linspace, arange, zeros, exp

import h5py
import networkx            as nx
import matplotlib.pyplot   as plt
import pygraphviz

#from matplotlib import rc
#rc('text', usetex=True)
#rc('font', size=14)
#rc('legend', fontsize=13)
#rc('text.latex', preamble=r'\usepackage{cmbright}')

def Directed_Graph(P):
    '''
        Dada la matriz de movilidad, esta función regresa el grafo con dirección y pesos
    '''
    # Create the graph with directional edges
    G = nx.DiGraph()
    
    for vi in range(N):
        for vj in range(N):
            
            if P[vi,vj] != 0: G.add_weighted_edges_from([(vi,vj,P[vi,vj])])
    
    return G
    
def Draw(G):

    pos = nx.circular_layout(G)

    plt.figure(figsize=(8,8))
    nx.draw(G,pos,node_color='none')

    nodes = nx.draw_networkx_nodes(G, pos,node_shape="o",node_size=700,alpha=1.0,node_color=[vi for vi in colors],with_labels=True)
    nodes.set_edgecolor('black')
    
    nodes_red = nx.draw_networkx_nodes(G,pos,node_size=750,node_shape="s",alpha=1.0,nodelist=[0],node_color='red',with_labels=True)
    nodes_red.set_edgecolor('black')

    nx.draw_networkx_labels(G,pos,labels = {i:i+1 for i in range(N) },font_color='white',font_size = 11 ) 

    plt.axis('equal')

    #if save == True: plt.savefig('Radial_Tree_Graph.png',bbox_inches='tight')
    plt.show()

# ====================================================================
#          ENTRADAS DEL PROGRAMA Y CONDICIONES INICIALES
# ====================================================================

# Condiciones iniciales
Initial_Conditions  = Models().Init_Conditions('uniform',Nh_range=[Nh_min, Nh_max])

# La matriz de movilidad
P           = Models().P(N,graph_name,pij_random)

# Arreglo con la distribuci\'on de valores beta (tasa de contacto )de cada parche
betas       = [random.uniform(beta_min,beta_max) for i in range(N)]

# Arreglo con la distribuci\'on de valores gamma (tasa de recuperaci\'on) de cada parche
gammas      = [random.uniform(gamma_min,gamma_max) for i in range(N)]

# Arreglo con los valores de la poblaci\'on total en cada parche
Ni          = [ Initial_Conditions[3*i] for i in range(N)  ]

G_directed  = Directed_Graph(P)

SOLVE       = True
ITERATE     = True
PLOT        = True

fname = graph_name + '_' + str(N)

# ================================================================================================================
#  RESOLVIENDO EL SISTEMA DE ECUACIONES CON EL CORRESPONDIENTE PROTOCOLO DE CONTROL Y GUARDANDO LA INFORMACI\'ON
# ================================================================================================================

if SOLVE == True:

    X   = odeint( Models().SIR, Initial_Conditions,time,args=(P,betas,gammas,Ni) )

    with h5py.File(fname+'.hdf5', 'w') as fh5:
        fh5.create_dataset("default",data=X)
        fh5.close()

# ================================================================================================================
#                                   RESOLVIENDO LAS ECUACIONES EN DIFERENCIAS
# ================================================================================================================

if ITERATE  == True:

    # Arreglo para cuando la condici\'on inicial es uno
    X_it_1 = zeros((N,int(px)))

    # Condici\'on inicial
    X_it_1[:,0] = 1.

    # Arreglo para cuando la condici\'on inicial es Ni
    X_it_Ni = zeros((N,int(px)))

    # Condici\'on inicial
    X_it_Ni[:,0] = Ni

    for step_i in range(1,int(px)):

        # Arreglo con el estado de los nodos en la iteraci\'on previa cuando la condici\'on inicial es uno
        X_prev_1  = X_it_1[:,step_i-1].reshape(1,N)[0]
    
        # Arreglo con el estado de los nodos en la iteraci\'on previa cuando la condici\'on inicial es S_i(0)
        X_prev_Ni = X_it_Ni[:,step_i-1].reshape(1,N)[0]
    
        # Ciclo sobre cada nodo
        for node_i in range(N):
            X_it_1[node_i][step_i]  = Ni[node_i] - Initial_Conditions[3*node_i]*exp(-Models().Theta_i(node_i,X_prev_1,P,betas,gammas,Ni))
            X_it_Ni[node_i][step_i] = Ni[node_i] - Initial_Conditions[3*node_i]*exp(-Models().Theta_i(node_i,X_prev_Ni,P,betas,gammas,Ni))

    # Plot format
    fontsize  = 15
    linewidth = 2
    labelsize = 15
    alpha     = 0.6
    formato   = 'png'

    file   = h5py.File(fname + '.hdf5','r')["default"]
    X      = file[:int(px/hpx)]

    #R_set  = [X[:,3*vi+2] for vi in range(N)]

    plt.figure(80,figsize=(15,7))
    
    cmap_1      = plt.cm.binary        #gist_rainbow   #Spectral #hsv
    colors_1    = cmap_1(linspace(0,1,N))
    
    cmap_2      = plt.cm.binary                 #gist_rainbow   #Spectral #hsv # gray
    colors_2    = cmap_2(linspace(0,1,N))

    #labeli = 'Patch %s' % str(ith_node +1)
    
    #plt.subplot(1,2,1)
    plt.plot(X_it_1[1], color = 'black',linewidth=linewidth,alpha = alpha) #,label=labeli
    plt.plot(X_it_1[10],color = 'gray',linewidth=linewidth,alpha = alpha)
    plt.plot(X_it_1[14],color = 'silver',linewidth=linewidth,alpha = alpha)
    
    plt.axhline(y=X[:,3*1+2][-1], linewidth=1.,color='red',linestyle='--',alpha=0.7)
    plt.axhline(y=X[:,3*10+2][-1],linewidth=1.,color='red',linestyle='--',alpha=0.7)
    plt.axhline(y=X[:,3*14+2][-1],linewidth=1.,color='red',linestyle='--',alpha=0.7)

    plt.ylabel("$I_{k}(t)$",fontsize=fontsize)
    plt.xlabel("$time$",fontsize=fontsize)
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    #plt.subplot(1,2,2)
    #plt.plot(X_it_Ni[1], color = colors_1[1],linewidth=linewidth,alpha = alpha)
    #plt.plot(X_it_Ni[10],color = colors_1[1],linewidth=linewidth,alpha = alpha)
    #plt.plot(X_it_Ni[14],color = colors_1[1],linewidth=linewidth,alpha = alpha)
    
    #plt.axhline(y=X[:,3*1+2][-1], linewidth=2.,color='red',linestyle='--')
    #plt.axhline(y=X[:,3*10+2][-1],linewidth=2.,color='red',linestyle='--')
    #plt.axhline(y=X[:,3*14+2][-1],linewidth=2.,color='red',linestyle='--')

    #plt.ylabel("$I_{k}(t)$",fontsize=fontsize)
    #plt.xlabel("$time$",fontsize=fontsize)
    #plt.tick_params(axis='y',labelsize=labelsize)
    #plt.tick_params(axis='x',labelsize=labelsize)

    #for ith_node in range(N):
        
    #    labeli = 'Patch %s' % str(ith_node +1)
        #plt.subplot(1,2,1)
    #    plt.plot(X_it_1[ith_node],color = colors_1[ith_node],linewidth=linewidth,alpha = alpha,label=labeli)
    #    plt.axhline(y=X[:,3*ith_node+2][-1],linewidth=2.,color='red',linestyle='--')

        #plt.subplot(1,2,2)
        #plt.plot(X_it_Ni[ith_node],color = colors_2[ith_node],linewidth=linewidth,alpha = alpha,label=labeli)
        #plt.axhline(y=X[:,3*ith_node+2][-1],linewidth=2.,color='red',linestyle='--')
    
    plt.show()

# ================================================================================================================
#                                    PLOT
# ================================================================================================================

if PLOT   == True:

    import numpy as np

    file   = h5py.File(fname + '.hdf5','r')["default"]
    X      = file[:int(px/hpx)]

    S_set  = [X[:,3*vi]   for vi in range(N)]
    I_set  = [X[:,3*vi+1] for vi in range(N)]
    R_set  = [X[:,3*vi+2] for vi in range(N)]

    Smean = zeros(int(px/hpx))
    for vi in range(int(px/hpx)):
        Smean[vi] = sum([S_set[vj][vi] for vj in range(N)])/N

    Imean = zeros(int(px/hpx))
    for vi in range(int(px/hpx)):
        Imean[vi] = sum([I_set[vj][vi] for vj in range(N)])/N

    Rmean = zeros(int(px/hpx))
    for vi in range(int(px/hpx)):
        Rmean[vi] = sum([R_set[vj][vi] for vj in range(N)])/N

    # Formato del grafico
    fontsize  = 15
    linewidth = 1.
    labelsize = 15
    alpha     = 0.6
    formato   = 'png'

    cmap       = plt.cm.binary                              #.gray #gist_rainbow #Spectral #hsv #jet
    colors     = cmap(linspace(0,1.0,N))

    # 
    #  NetworkX graph_name
    #

    fig =plt.figure(10,figsize=(20,20),dpi=80)   #, facecolor='w', edgecolor='k')
    #fig.set_size_inches(10.0, 20.0, 10)

    plt.title('a)  Networked popution',fontsize=26)
    #plt.text(-1.1, 1.1, r'a)', horizontalalignment='left',verticalalignment='top',fontsize=36)


    pos = nx.circular_layout(G_directed)

    #edge_labels=dict([((u,v,),d['weight']) for u,v,d in G_directed.edges(data=True)])
    #nx.draw_networkx_edge_labels(G_directed,pos,edge_labels=edge_labels)

    nx.draw(G_directed,pos,node_color='none')

    nx.draw_networkx_edges(G_directed, pos, width=0.5,arrowsize=20)

    nodes = nx.draw_networkx_nodes(G_directed, pos,node_shape="o",node_size=600,edgecolors='black',alpha=1.0,node_color='gray',with_labels=True,linewidths=1.5) #[colors[vi] for vi in range(N)]
    nodes.set_edgecolor('black')
    
    nodes_red = nx.draw_networkx_nodes(G_directed,pos,node_size=650,node_shape="s",alpha=1.0,nodelist=[0],node_color='gray',with_labels=True,linewidths=1.5)
    nodes_red.set_edgecolor('black')

    nx.draw_networkx_labels(G_directed,pos,labels = {i:i+1 for i in range(N) },font_color='white',font_size = 18 ) 

    plt.axis('equal')

    #plt.savefig('NetworkX.png')

    #plt.show()

    #
    #  Matrix P
    #

    fig =plt.figure(20,figsize=(20,20),dpi=80)

    plt.title('b) Residence-time Matrix',fontsize=26)

    im = plt.imshow(P, interpolation='nearest') #cmap=plt.cm.gray
    for i in range(N):
        for j in range(N):
            text = plt.text(j, i, round(P[i, j],1),ha="center", va="center", color="w",fontsize=14)

    plt.xticks( [i for i in range(N)],labels = [str(i+1) for i in range(N)] )
    plt.yticks( [i for i in range(N)],labels = [str(i+1) for i in range(N)] )

    plt.tick_params(axis='both', which='major', labelsize=18)

    plt.xlabel(r'Patch label',fontsize=fontsize+10)
    plt.ylabel(r'Patch label',fontsize=fontsize+10)

    #plt.savefig('Matrix.png')

    #plt.show()

    #
    #  Infected
    #

    fig =plt.figure(30,figsize=(20,20),dpi=80)
    plt.subplot(1,3,1)
    plt.title('a) Number of susceptibles individuals in each patch',fontsize=15)
    #plt.subplot(1,3,3)
    for vi in range(N):
        plt.plot(time, X[:,3*vi], color=colors[vi],linewidth=2.,alpha=0.5,label='Patch ' + str(vi+1))
    
    plt.plot(time,Smean,'--',linewidth=5.,color='red')

    plt.xlabel(r'time (days)',fontsize=fontsize)
    #plt.ylabel(r'Number of susceptibles',fontsize=fontsize+10)
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.xlim([0,15])

    #plt.savefig('susceptibles.png')

    #fig =plt.figure(40,figsize=(20,20),dpi=80)
    plt.subplot(1,3,2)
    plt.title('b) Number of infected individuals in each patch',fontsize=15)
    #plt.subplot(1,3,3)
    for vi in range(N):
        plt.plot(time, X[:,3*vi+1], color=colors[vi],linewidth=2.,alpha=0.5,label='Patch ' + str(vi+1))
    
    plt.plot(time,Imean,'--',linewidth=5.,color='red')

    plt.xlabel(r'time (days)',fontsize=fontsize)
    #plt.ylabel(r'Number of infected',fontsize=fontsize+10)
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.xlim([0,15])
    #plt.ylim([0:])

    #plt.legend (bbox_to_anchor=(0.98,0.98),loc=1,borderaxespad=0.)

    #plt.savefig('infected.png')

    #fig =plt.figure(50,figsize=(20,20),dpi=80)
    plt.subplot(1,3,3)
    plt.title('c) Number of recovered individuals in each patch',fontsize=15)
    #plt.subplot(1,3,3)
    for vi in range(N):
        plt.plot(time, X[:,3*vi+2], color=colors[vi],linewidth=2.,alpha=0.5,label='Patch ' + str(vi+1))
    
    plt.plot(time,Rmean,'--',linewidth=5.,color='red')

    plt.xlabel(r'time (days)',fontsize=fontsize)
    #plt.ylabel(r'Number of recovered ',fontsize=fontsize+10)
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.xlim([0,15])
    #plt.ylim([0:])

    #plt.legend (bbox_to_anchor=(0.98,0.98),loc=1,borderaxespad=0.)
    #plt.savefig('recovered.png')

    plt.show()
