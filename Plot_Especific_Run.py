# -*- coding: utf-8 -*-
"""
@name  : Plot_Especific_Run.py
@author: ANDRES
@date:   JULY 28, 2020
@update: 

PROGRAMA PARA REALIZAR LAS GRAFICAS DEL ARTICULO PARA LOS SIGUIENTES EJEMPLOS SELECCIONADOS

GRAFO      RUN PARA ABSOLUTO     RUN PARA RELATIVO
ER               33                     99
SWN              9                       4
BA               100                    46

"""

import h5py, pickle, sys
from Settings             import *
from numpy                import linspace, argmax
import matplotlib.pyplot   as plt

from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter


def num_format(number):
    return ("{:,}".format(round(number,2)))

# =================================================================================================
#                               SETTINGS
# =================================================================================================

Folders_Names = ['NWS_Absoluto_RUN_9','NWS_Relativo_RUN_4','BA_Absoluto_RUN_100','BA_Relativo_RUN_46','ER_Absoluto_RUN_33','ER_Relativo_RUN_99']

#Folders_Names = ['BA_Absoluto_RUN_100']

# Formato del grafico
fontsize  = 20
linewidth = 1.
labelsize = 15
alpha     = 0.6
formato   = 'png'

cmap      = plt.cm.binary                  #.gray #gist_rainbow #Spectral #hsv #jet
colors    = cmap(linspace(0,1.0,N))

# Font size for the extra text with the maximum of I
font_I_data = 20

for cont,File_Name in enumerate(Folders_Names):

    # =================================================================================================
    #                               OPEN H5 FILES
    # =================================================================================================

    # The run number of the selected patch.
    run_number = File_Name.split('_')[-1]

    # Abriendo el archivo pkl con la información de cuales fueron los nodos controlados
    infile         = open( File_Name + '/M_SELECTED_PATCH_RUN_' + run_number + '.pkl','rb')
    new_dict       = pickle.load(infile)

    # El conjunto de parches seleccionados con el índice propuesto
    ctr_patch      = new_dict[0:1]
    selected_patch = ctr_patch[0]

    # El conjunto de parches seleccionados aleatoriamente  --> En esta versión considero que solo un parche es seleccionado
    selected_rand  = new_dict[1]

    # Abriendo el archivo hdf5 para la simulacion sin control
    file_no_Ctr    = h5py.File(File_Name + '/Without_Control_Run_' + run_number +'.hdf5','r')["default"]
    X_no_Ctr       = file_no_Ctr[:int(px/hpx)]

    # Abriendo el archivo hdf5 para la simulacion controlando el parche con alguno de los índices (relativo o absoluto) - Protocolo 1
    file_Ctr_Max_I = h5py.File(File_Name + '/Max_I_Run_' + run_number +'.hdf5','r')["default"]
    X_Ctr_Max_I    = file_Ctr_Max_I[:int(px/hpx)]

    # Abriendo el archivo hdf5 para la simulacion controlando aleatoriamente el parche.
    file_Ctr_Random = h5py.File(File_Name + '/Random_Patches_Run_' +  run_number +'.hdf5','r')["default"]
    X_Ctr_Random    = file_Ctr_Random[:int(px/hpx)]

    # =================================================================================================
    #                               PLOT
    # =================================================================================================

    fig = plt.figure(101+cont,figsize=(20,10))
    plt.subplots_adjust(bottom=0.1, right=0.98, left=0.09, top=0.9,wspace=0.2, hspace=0.3)
    #fig.subplots_adjust(wspace=0.2, hspace=0.3)
    #ax  = fig.add_subplot(111, autoscale_on=False, xlim=(-1, 5), ylim=(-3, 5))

    # =========================================================================
    #                   SUB-PLOT 1
    # =========================================================================

    plt.subplot(3,3,1)
    plt.title('Susceptibles',fontsize=20)

    plt.annotate('(a)', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.95, 0.95), textcoords='axes fraction',fontsize=20,
    horizontalalignment='right', verticalalignment='top')

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
    plt.yticks([200000,450000,600000,800000,1000000], [0.2,0.4,0.6,0.8,1],rotation='horizontal') 

    plt.annotate(r'$\times 10^{6}$', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.09, 1.1), textcoords='axes fraction',fontsize=15,
    horizontalalignment='right', verticalalignment='top')
     
    #plt.margins(0.2)
    plt.xlim([0,20])
    plt.ylim([0,1000000])

    # =========================================================================
    #                   SUB-PLOT 2
    # =========================================================================

    plt.subplot(3,3,2)
    plt.title('Infected',fontsize=20)

    plt.annotate('(b)', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.1, 0.95), textcoords='axes fraction',fontsize=20,
    horizontalalignment='right', verticalalignment='top')

    for vi in range(N):
        plt.plot(time, X_no_Ctr[:,3*vi+1], color=colors[vi],linewidth=2.,alpha=0.5)
    
    # To plot the selected patch with color red
    for select_i in ctr_patch:
        plt.plot(time, X_no_Ctr[:,3*select_i+1], color = 'red',linewidth=2.5,alpha=1.0)

    # Tiempo para llegar al máximo de los infectados
    time_max_select_i = argmax(X_no_Ctr[:,3*selected_patch+1])*hpx

    plt.scatter(time_max_select_i,max(X_no_Ctr[:,3*selected_patch+1]),marker='o',color='red',s=30)
    #plt.scatter(time_max_select_i,max(X_no_Ctr[:,3*selected_rand+1]),marker='o',color='green',s=30)
    plt.scatter(4.931,max(X_no_Ctr[:,3*selected_rand+1]),marker='o',color='green',s=30)

    plt.annotate( num_format( max(X_no_Ctr[:,3*selected_patch+1]) ), xy=(time_max_select_i, max(X_no_Ctr[:,3*selected_patch+1])),  xycoords='data',
    xytext=(0.9, 0.95), textcoords='axes fraction',fontsize=font_I_data,
    horizontalalignment='right', verticalalignment='top',
    arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})
    
    # To plot the selected patch selected randomly with color green
    plt.plot(time, X_no_Ctr[:,3*selected_rand+1],color = 'green',  linewidth=2.5,alpha=1.0)

    plt.annotate( num_format( max(X_no_Ctr[:,3*selected_rand+1]) ), xy=(time_max_select_i, max(X_no_Ctr[:,3*selected_rand+1])),  xycoords='data',
    xytext=(0.36, 0.75), textcoords='axes fraction',fontsize=font_I_data,
    horizontalalignment='right', verticalalignment='top',
    arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)
    plt.xlim([0,20])
    plt.ylim([0,400000])

    #plt.yticks([100000,200000,300000,400000], [0.1,0.2,0.3,0.4],rotation='horizontal')
    plt.yticks([100000,200000,300000,400000,500000], [0.1,0.2,0.3,0.4,0.5],rotation='horizontal') 

    plt.annotate(r'$\times 10^{6}$', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.09, 1.1), textcoords='axes fraction',fontsize=15,
    horizontalalignment='right', verticalalignment='top')

    # =========================================================================
    #                   SUB-PLOT 3
    # =========================================================================

    plt.subplot(3,3,3)
    plt.title('Recovered',fontsize=20)

    plt.annotate('(c)', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.1, 0.95), textcoords='axes fraction',fontsize=20,
    horizontalalignment='right', verticalalignment='top')

    for vi in range(N):
        plt.plot(time, X_no_Ctr[:,3*vi+2], color=colors[vi],linewidth=2.,alpha=0.5)
    
    # To plot the selected patch with color red
    for select_i in ctr_patch:
        plt.plot(time, X_no_Ctr[:,3*select_i+2], color = 'red',linewidth=2.5,alpha=1.0)
    
    plt.annotate( num_format( max(X_no_Ctr[:,3*selected_patch+2]) ), xy=(17.5, max(X_no_Ctr[:,3*selected_patch+2])),  xycoords='data',
    xytext=(0.5, 0.95), textcoords='axes fraction',fontsize=font_I_data,
    horizontalalignment='right', verticalalignment='top',
    arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

    # To plot the selected patch selected randomly with color green
    plt.plot(time, X_no_Ctr[:,3*selected_rand+2],color = 'green',  linewidth=2.5,alpha=1.0)

    plt.annotate( num_format( max(X_no_Ctr[:,3*selected_rand+2]) ), xy=(17.5, max(X_no_Ctr[:,3*selected_rand+2])),  xycoords='data',
    xytext=(0.38, 0.5), textcoords='axes fraction',fontsize=font_I_data,
    horizontalalignment='right', verticalalignment='top',
    arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.yticks([200000,400000,600000,800000,1000000], [0.2,0.4,0.6,0.8,1],rotation='horizontal')

    plt.annotate(r'$\times 10^{6}$', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.09, 1.1), textcoords='axes fraction',fontsize=15,
    horizontalalignment='right', verticalalignment='top')

    plt.xlim([0,20])
    plt.ylim([0,1000000])

    # =========================================================================
    #                   SUB-PLOT 4
    # =========================================================================

    plt.subplot(3,3,4)

    plt.annotate('(d)', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.95, 0.95), textcoords='axes fraction',fontsize=20,
    horizontalalignment='right', verticalalignment='top')
    
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

    plt.yticks([200000,400000,600000,800000,1000000], [0.2,0.4,0.6,0.8,1],rotation='horizontal') 

    plt.annotate(r'$\times 10^{6}$', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.09, 1.1), textcoords='axes fraction',fontsize=15,
    horizontalalignment='right', verticalalignment='top')
     
    #plt.margins(0.2)
    plt.xlim([0,20])
    plt.ylim([0,1000000])

    # =========================================================================
    #                   SUB-PLOT 5
    # =========================================================================

    plt.subplot(3,3,5)

    plt.annotate('(e)', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.1, 0.95), textcoords='axes fraction',fontsize=20,
    horizontalalignment='right', verticalalignment='top')
    
    for vi in range(N):
        plt.plot(time, X_Ctr_Max_I[:,3*vi+1], color=colors[vi],linewidth=2.,alpha=0.5)
    
    # To plot the selected patch with color red
    for select_i in ctr_patch:
        plt.plot(time, X_Ctr_Max_I[:,3*select_i+1], color = 'red',linewidth=2.5,alpha=1.0)

    # Tiempo para llegar al máximo de los infectados
    time_max_select_max = argmax(X_Ctr_Max_I[:,3*selected_patch+1])*hpx

    plt.scatter(time_max_select_max,max(X_Ctr_Max_I[:,3*selected_patch+1]),marker='o',color='red',s=30)
    #plt.scatter(time_max_select_max,max(X_Ctr_Max_I[:,3*selected_rand+1]),marker='o',color='green',s=30)
    plt.scatter(5.074,max(X_Ctr_Max_I[:,3*selected_rand+1]),marker='o',color='green',s=30)

    plt.annotate( num_format( max(X_Ctr_Max_I[:,3*selected_patch+1]) ), xy=(time_max_select_max, max(X_Ctr_Max_I[:,3*selected_patch+1])),  xycoords='data',
    xytext=(0.9, 0.95), textcoords='axes fraction',fontsize=font_I_data,
    horizontalalignment='right', verticalalignment='top',
    arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

    # To plot the selected patch selected randomly with color green
    plt.plot(time, X_Ctr_Max_I[:,3*selected_rand+1],color = 'green',  linewidth=2.5,alpha=1.0)

    plt.annotate( num_format( max(X_Ctr_Max_I[:,3*selected_rand+1]) ), xy=(time_max_select_max, max(X_Ctr_Max_I[:,3*selected_rand+1])),  xycoords='data',
    xytext=(0.36, 0.75), textcoords='axes fraction',fontsize=font_I_data,
    horizontalalignment='right', verticalalignment='top',
    arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    #plt.yticks([100000,200000,300000,400000], [0.1,0.2,0.3,0.4],rotation='horizontal') 
    plt.yticks([100000,200000,300000,400000,500000], [0.1,0.2,0.3,0.4,0.5],rotation='horizontal')

    plt.annotate(r'$\times 10^{6}$', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.09, 1.1), textcoords='axes fraction',fontsize=15,
    horizontalalignment='right', verticalalignment='top')

    plt.xlim([0,20])
    plt.ylim([0,400000])

    # =========================================================================
    #                   SUB-PLOT 6
    # =========================================================================

    plt.subplot(3,3,6)

    plt.annotate('(f)', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.1, 0.95), textcoords='axes fraction',fontsize=20,
    horizontalalignment='right', verticalalignment='top')
    
    for vi in range(N):
        plt.plot(time, X_Ctr_Max_I[:,3*vi+2], color=colors[vi],linewidth=2.,alpha=0.5)
    
    # To plot the selected patch with color red
    for select_i in ctr_patch:
        plt.plot(time, X_Ctr_Max_I[:,3*select_i+2], color = 'red',linewidth=2.5,alpha=1.0)
    
    plt.annotate( num_format( max(X_Ctr_Max_I[:,3*selected_patch+2]) ), xy=(17.5, max(X_Ctr_Max_I[:,3*selected_patch+2])),  xycoords='data',
    xytext=(0.5, 0.95), textcoords='axes fraction',fontsize=font_I_data,
    horizontalalignment='right', verticalalignment='top',
    arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})
    
    # To plot the selected patch selected randomly with color green
    plt.plot(time, X_Ctr_Max_I[:,3*selected_rand+2],color = 'green',  linewidth=2.5,alpha=1.0)

    plt.annotate( num_format( max(X_Ctr_Max_I[:,3*selected_rand+2]) ), xy=(17.5, max(X_Ctr_Max_I[:,3*selected_rand+2])),  xycoords='data',
    xytext=(0.4, 0.5), textcoords='axes fraction',fontsize=font_I_data,
    horizontalalignment='right', verticalalignment='top',
    arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.yticks([200000,400000,600000,800000,1000000], [0.2,0.4,0.6,0.8,1],rotation='horizontal')

    plt.annotate(r'$\times 10^{6}$', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.09, 1.1), textcoords='axes fraction',fontsize=15,
    horizontalalignment='right', verticalalignment='top')

    plt.xlim([0,20])
    plt.ylim([0,1000000])
    
    # =========================================================================
    #                   SUB-PLOT 7
    # =========================================================================

    plt.subplot(3,3,7)

    plt.annotate('(g)', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.95, 0.95), textcoords='axes fraction',fontsize=20,
    horizontalalignment='right', verticalalignment='top')

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

    plt.yticks([200000,400000,600000,800000,1000000], [0.2,0.4,0.6,0.8,1],rotation='horizontal') 

    plt.annotate(r'$\times 10^{6}$', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.09, 1.1), textcoords='axes fraction',fontsize=15,
    horizontalalignment='right', verticalalignment='top')
     
    #plt.margins(0.2)
    plt.xlim([0,20])
    plt.ylim([0,1000000])

    # =========================================================================
    #                   SUB-PLOT 8
    # =========================================================================

    plt.subplot(3,3,8)

    plt.annotate('(h)', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.1, 0.95), textcoords='axes fraction',fontsize=20,
    horizontalalignment='right', verticalalignment='top')

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
    xytext=(0.9, 0.95), textcoords='axes fraction',fontsize=font_I_data,
    horizontalalignment='right', verticalalignment='top',
    arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})
    
    # To plot the selected patch selected randomly with color green
    plt.plot(time, X_Ctr_Random[:,3*selected_rand+1],color = 'green',  linewidth=2.5,alpha=1.0)

    plt.annotate( num_format( max(X_Ctr_Random[:,3*selected_rand+1]) ), xy=(time_max_select_rand, max(X_Ctr_Random[:,3*selected_rand+1])),  xycoords='data',
    xytext=(0.36, 0.75), textcoords='axes fraction',fontsize=font_I_data,
    horizontalalignment='right', verticalalignment='top',
    arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

    plt.xlabel(r'time (days)',fontsize=fontsize)
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)
    plt.xlim([0,20])
    plt.ylim([0,400000])

    #plt.yticks([100000,200000,300000,400000], [0.1,0.2,0.3,0.4],rotation='horizontal') 
    plt.yticks([100000,200000,300000,400000,500000], [0.1,0.2,0.3,0.4,0.5],rotation='horizontal')

    plt.annotate(r'$\times 10^{6}$', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.09, 1.1), textcoords='axes fraction',fontsize=15,
    horizontalalignment='right', verticalalignment='top')

    # =========================================================================
    #                   SUB-PLOT 9
    # =========================================================================

    plt.subplot(3,3,9)

    plt.annotate('(i)', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.1, 0.95), textcoords='axes fraction',fontsize=20,
    horizontalalignment='right', verticalalignment='top')

    for vi in range(N):
        plt.plot(time, X_Ctr_Random[:,3*vi+2], color=colors[vi],linewidth=2.,alpha=0.5)
    
    # To plot the selected patch with color red
    for select_i in ctr_patch:
        plt.plot(time, X_Ctr_Random[:,3*select_i+2], color = 'red',linewidth=2.5,alpha=1.0)
    
    plt.annotate( num_format( max(X_Ctr_Random[:,3*selected_patch+2]) ), xy=(17.5, max(X_Ctr_Random[:,3*selected_patch+2])),  xycoords='data',
    xytext=(0.5, 0.95), textcoords='axes fraction',fontsize=font_I_data,
    horizontalalignment='right', verticalalignment='top',
    arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

    # To plot the selected patch selected randomly with color green
    plt.plot(time, X_Ctr_Random[:,3*selected_rand+2],color = 'green',  linewidth=2.5,alpha=1.0)

    plt.annotate( num_format( max(X_Ctr_Random[:,3*selected_rand+2]) ), xy=(17.5, max(X_Ctr_Random[:,3*selected_rand+2])),  xycoords='data',
    xytext=(0.38, 0.5), textcoords='axes fraction',fontsize=font_I_data,
    horizontalalignment='right', verticalalignment='top',
    arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"})

    plt.xlabel(r'time (days)',fontsize=fontsize)
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.yticks([200000,400000,600000,800000,1000000], [0.2,0.4,0.6,0.8,1],rotation='horizontal')

    plt.annotate(r'$\times 10^{6}$', xy=(1.1, 2.1),  xycoords='data',
    xytext=(0.09, 1.1), textcoords='axes fraction',fontsize=15,
    horizontalalignment='right', verticalalignment='top')

    plt.xlim([0,20])
    plt.ylim([0,1000000])

    plt.savefig(File_Name + '.png')

    del time_max_select_i,time_max_select_max, time_max_select_rand

    #plt.show()
