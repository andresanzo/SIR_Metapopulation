"""
@name  : Run_SIMS.py
@author: ANDRES
@date:   MARZO 18, 2020
@update: ABRIL 7,  2020

EL PROGRAMA REALIZA LOS SIGUIENTES PASOS POR CADA SIMULACI\'ON

STEP 0:
        -- GENERA LA RED DE MOVILIDAD CON LA ESTRUTURA DEFINIDA EN EL ARCHIVO Settings.py
        -- GENERA LAS CONDICIONES INICIALES
        -- GENERA LA DISTRIBUCI\'ON DE BETAS
        -- TODO LO ANTERIOR LO REALIZA EJECUTANDO EL ARCHIVO Generate_Initial_Conditions.py  <Numero de Run>
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
       -- REALIZAR LA SIMULACI\'O CONTROLANDO M PARCHES SELECCIONADOS DE FORMA ALEATORIA
       
STEP 6:
     -- REALIZAR SIMULACI\'ON SIN CONTROL

STEP 7:
      -- GUARDAR LA INFORMACI\'ON EN LA CARPETA RUN_<run_i>
      
    
 ENTRADAS DEL SISTEMA
    --

 LA SALIDA DEL SISTEMA
    -- UNA CARPERTA POR CADA CORRIDA CON LAS SERIES DE TIEMPO DE LAS SIMULACIONES
    
Para compilar este programa desde una consala teclear  python3 Run_SIMS.py <Numero de parches a controlar M>

"""

import pickle,csv,sys
from numpy         import random
from subprocess    import call
from Settings      import *
from SIR_Functions import Models
from SIR_Functions import Save_Info
from tqdm          import tqdm
import networkx    as     nx
from heapq import nlargest   # Para hacer la selecci\'on de los M parches con el max I_max

# El Número de parches a controlar. 
M = int(sys.argv[1])

# Ciclo sobre cada corrida
for run_i in tqdm(range(N_runs_ensamble_min,N_runs_ensamble_max)):


        # ========================================================================================
        #                                       STEP 0
        # ========================================================================================

        # Genera el archivo 'INITIAL_CONDITIONS_RUN_<Run_Number>.pkl con la informaci\'on de la matriz de movilidad, condiciones iniciales y betas
        call(['python3','Generate_Initial_Conditions.py',str(run_i)])
    
        # Abriendo el archivo pkl
        infile   = open('INITIAL_CONDITIONS_RUN_' + str(run_i) +'.pkl','rb')
        new_dict = pickle.load(infile)
        
        # Condiciones iniciales
        Initial_Conditions = new_dict[0]

        # La matriz de movilidad
        P                  = new_dict[1]

        # Los valores beta para cada parche
        betas              = new_dict[2]
    
        # Los valores beta para cada parche
        gammas             = new_dict[3]

        # Los valores de la poblaci\'on total en cada parche
        Ni                 =  new_dict[4]
        
        # La matriz de movilidad
        G                  = new_dict[5]
        
        # ========================================================================================================================
        #                                       STEP 1
        # ========================================================================================================================
    
        # Resuleve el sistema en diferencias y regresa la lista de I_total de cada parche
        I_total   = Models().Controlled_Patch(Initial_Conditions,select_protocol,P,betas,gammas,Ni)
        
        #print('----->',run_i,'<-------')
        #print(I_total)
        #print('------------------------')
    
        # El parche con el maximo porcentaje de infectados
        #ctr_patch = I_total.index(max(I_total))
        
        # Selecciona los M parches que de acuerdo al \'indice, tienen el max num. de infectados. El arreglo contiene índice de los parches
        ctr_patches = sorted(range(len(I_total)), key = lambda sub: I_total[sub])[-M:]
        
        # ========================================================================================================================
        #                                       STEP 6
        # ========================================================================================================================
        
        # Realiza la simulaci\'on sin control
        call(['python3','SIR_Controlled_M_Patches.py','[-1]',str(run_i),'Without_Control_Run_'+str(run_i),'Sin_control'])
        
        # ========================================================================================================================
        #                                       STEP 2
        # ========================================================================================================================
    
        call(['python3','SIR_Controlled_M_Patches.py',str(ctr_patches),str(run_i),'Max_I_Run_'+str(run_i),'Control_Beta'])
    
        # ========================================================================================================================
        #                                       STEP 3
        # ========================================================================================================================
    
        #call(['python3','SIR_Controlled_M_Patches.py',str(ctr_patches),str(run_i),'Pij_in_Run_'+str(run_i),'Control_Pij_in'])
        
        # ========================================================================================================================
        #                                       STEP 4
        # ========================================================================================================================
    
        #call(['python3','SIR_Controlled_M_Patches.py',str(ctr_patches),str(run_i),'Pij_out_Run_'+str(run_i),'Control_Pij_out'])

        # ========================================================================================================================
        #                                       STEP 5
        # ========================================================================================================================
    
        # Selecciona de forma aleatoria M parches ( con M definido en Settings.py )
        
        # Las etiquetas de los parches sin los M seleccionados de acuerdo al índice
        Labels = [i for i in range(N) if i not in ctr_patches]
        
        random_patches = []
        for vi in range(M):
            
            random_select = random.choice(Labels)
            
            while  random_select in random_patches:
                random_select   = random.choice(Labels)
                
            random_patches.append(random_select)
        
        outfile = open('M_SELECTED_PATCH_RUN_' + str(run_i) + '.pkl','wb')
        pickle.dump(ctr_patches + random_patches,outfile)
        
        call(['python3','SIR_Controlled_M_Patches.py',str(random_patches),str(run_i),'Random_Patches_Run_'+str(run_i),'Control_Beta'])
    
        
        # ========================================================================================================================
        #                                       STEP 7
        # =========================================================================================================================
        
        # Crea la carpeta RUN_<run_i>
        Save_Info().MakeDir(graph_name + '_M_' + str(M) + '/RUN_'+str(run_i))
    
        # Guarda el archivo pkl generado en la carpeta RUN_<run_i>
        Save_Info().SaveFiles(graph_name + '_M_' + str(M) + '/RUN_'+str(run_i),make_dir=False,file_format='*.pkl')
    
        # Guarda los archivos hdf5 generado en la carpeta RUN_<run_i>
        Save_Info().SaveFiles(graph_name + '_M_' + str(M) + '/RUN_'+str(run_i),make_dir=False,file_format='*.hdf5')
        
        infile.close()
        outfile.close()

