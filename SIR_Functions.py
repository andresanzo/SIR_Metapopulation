# -*- coding: utf-8 -*-
# DATE: FEBRERO 22, 2020.

#import random
from __future__ import division

import shutil,os,glob
import networkx             as     nx
import matplotlib.pyplot    as     plt

from   numpy                import zeros,random,array,full_like,transpose,exp
from   scipy                import arange
from   Settings             import *        

class Models:
    
    def __init__(self):
        return
        
    def StringArray_to_FloatArray(self,A):
        ''' Program to transform an array in string format to an array of float items
            A  ---> String of the form [ [a,b],[c,d],...,[y,z]]
        '''
        a = A.replace('[','').replace(']','').split(',')

        return array([int(ai) for ai in a])
        
    def Init_Conditions(self,Nh_range=[10000,1000000]):
        
        """
        Function to generate randomly the initial condition for each patch
        """
        
        init_cond = zeros(3*N)

        for i in range(N):
            init_cond[3*i]   = random.uniform(Nh_range[0],Nh_range[1])  # ------------> Number of host Nh
        
        # The node with label zero is the selected node where a single infected is introduced
        init_cond[1] = 1.0

        # The distribution of population size in each patch
        Nhi = [init_cond[3*k] for k in range(N)]

        return init_cond, Nhi
        
    def P(self,N,graph_name,pij_random):

        """
            Program to generate the dwell-time matrix P

            N           ---> Float  ---> Number of patches
            graph_name  ---> String ---> Option of network topology
            pij_random  ---> String ---> Option of the type of pij random selection

        """
        
        if graph_name == 'Fully': G = nx.complete_graph(N)
        if graph_name == 'ER'   : G = nx.erdos_renyi_graph(N,0.6)
        if graph_name == 'NWS'  : G = nx.newman_watts_strogatz_graph(N,4,0.3)
        if graph_name == 'BA'   : G = nx.barabasi_albert_graph(N,2)

        # CONSTRUCT THE P MATRIX
        if pij_random == 'uniform':
            P = zeros((N,N))
            for v in G.edges(): 
                P[v[0],v[1]] = random.uniform(0.1,0.9)
                P[v[1],v[0]] = random.uniform(0.1,0.9)
    
            # THE DIAGONAL ENTRIES OF P 
            for i in range(N): P[i,i] =  random.uniform(0.1,0.9)
        
            # NORMALIZED P
            P_normed = array([[P[i,j]/float(sum(P[i])) for j in range(N)] for i in range(N) ])
    
        if pij_random == 'power':
            P = zeros((N,N))
            for v in G.edges(): 
                P[v[0],v[1]] = random.power(1.7,1) 
                P[v[1],v[0]] = random.power(1.7,1) 
    
            # THE DIAGONAL ENTRIES OF P 
            for i in range(N): P[i,i] =  random.power(1.7,1) 

            # NORMALIZED P
            P_normed = array([[P[i,j]/float(sum(P[i])) for j in range(N)] for i in range(N) ])

        return P_normed,G
        
    def K_in(self,P):
        '''
        Regresa un arrelgo con la lista del grado de nodo para los enlaces hacia el interior del parche k_in.
        
        Entrada:
            P --- Arreglo: Matriz de movildad
        Salida
           Arreglo con los valores k_in de cada nodo
        '''
        return [sum(P[:,node_i]) - P[node_i,node_i]  for node_i in range(N) ]
        
    def K_out(self,P):
        '''
        Regresa un arrelgo con la lista del grado de nodo para los enlaces hacia el exterior del parche k_out.
    
        Entrada:
            P --- Arreglo: Matriz de movildad
        Salida
        Arreglo con los valores k_out de cada nodo
        '''
        return [sum(P[node_i,:]) - P[node_i,node_i]  for node_i in range(N)]
        
    def weighted_choice(self,weights):
    
        """
        The following is a simple function to implement weighted random selection in Python. Given a list of weights,   
        it returns an index randomly, according to these weights. For example, given [2, 3, 5] it returns 0 (the index  
        of the first element) with probability 0.2, 1 with probability 0.3 and 2 with probability 0.5. The weights  
        need    
        not sum up to anything in particular, and can actually be arbitrary Python floating point numbers.
        http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/#id1
        """

        totals = []
        running_total = 0

        for w in weights:
            running_total += w
            totals.append(running_total)

        rnd = random.random() * running_total
        for i, total in enumerate(totals):
            if rnd < total:
                return i
                
    def SIR(self,X,t,P,betas,gammas,Ni):

        # Main Loop to compute the nodes and links weights dynamical state
        odes = []

        for i in range(N):
        
            Sum = sum([ self.alpha_ij(i,j,P,betas,Ni)*X[3*j+1] for j in range(N)] )
            
            odes.extend([
                          - X[3*i]*Sum ,
                            X[3*i]*Sum  - gammas[i]*X[3*i + 1] ,
                            gammas[i]*X[3*i + 1]
                        ])
                
        return odes
                
    #def W_h(self,X,P,j):
    #    return sum([(X[3*ki]+X[3*ki+1]+X[3*ki+2])*P[ki,j] for ki in range(N)])
        
    def W_h(self,Ni,P,i):
        ''' Poblaci\'on efectiva en el parche i
            Ni   -- Arreglo: los valores de la poblaci'on total en cada parche
            P    -- Arreglo: matriz de movilidad
            i    -- int: \'indice del parche '''
        return sum([Ni[ki]*P[ki,i] for ki in range(N)])
        
    def alpha_ij(self,i,j,P,betas,Ni):
        ''' Calcula el t\'ermino alpha_ij del modelo para el par de parches con etiquetas i y j.
            i,j   -- float: equiquetas de los parches
            P     -- Arreglo: matriz de movilidad
            betas -- Arreglo: valores de beta para cada parche '''
        return sum([ (betas[k]*P[i,k]*P[j,k])/(self.W_h(Ni,P,j))  for k in range(N)  ])
        
    def Theta_i(self,i,X,P,betas,gammas,Ni):
        ''' Calcula el \'indice Theta para el parche con etiqueta i
            i     -- float: equiqueta del parche
            X     -- Arreglo: Estados (Infectados) de los parches
            betas -- Arreglo: valores de beta para cada parche
            Ni    -- Arreglo: Distribuci\'on de la densididad poblacional de cada parche '''
        return sum( [ (self.alpha_ij(i,j,P,betas,Ni)*X[j])/(gammas[j]) for j in range(N) ]  )
        
    def Iterate_Difference_System(self,Init_Conds,P,betas,gammas,Ni):
        ''' Itera el sistema en diferencias X_{n+1} = Ni - S_i(0)exp{-Theta_i}
            X     -- Arreglo: Estados (Infectados) de los parches
            betas -- Arreglo: valores de beta para cada parche
            Ni    -- Arreglo: Distribuci\'on de la densididad poblacional de cada parche
        '''
        # N\'umero de iteraciones
        num_steps = 1000
        
        # Arreglo para guardar las series de tiempo (N es el n\'umero de nodos definido en Settings.py
        X_it = zeros((N,num_steps))
        
        # Condici\'on inicial
        X_it[:,0] = Ni
        
        # Iterando
        for step_i in range(1,num_steps):
            
            X_prev = X_it[:,step_i-1].reshape(1,N)[0]
            
            # Ciclo sobre cada nodo
            for node_i in range(N):
                X_it[node_i][step_i] = Ni[node_i] - Init_Conds[3*node_i]*exp(-self.Theta_i(node_i,X_prev,P,betas,gammas,Ni))
        
        return X_it
        
    def Controlled_Patch(self,Init_Conds,selection_protocol,P,betas,gammas,Ni):
        ''' Dada la matriz de movilidad, la distribuci\'on de betas y densidad poblacional, esta funci\'on regresa el parche con el
            maximo porcentaje de infectados de acuerdo con el sistema en diferencias.
            seletion_protol (string) tienen dos opciones: 'relative' para sleccionar el parche con el max I_i/Ni y 'absolute' para seleccionar max I_i
        '''
        
        # Resuelve el sistema
        X = self.Iterate_Difference_System(Init_Conds,P,betas,gammas,Ni)
        
        if selection_protocol == 'relative':
            # El tamaño final de cada parche dividido entre Ni
            I_total = [ X[:,-1][node_i]/Ni[node_i] for node_i in range(N)]
            
        if selection_protocol == 'absolute':
            # El tamaño final de cada parche dividido entre Ni
            I_total = [ X[:,-1][node_i] for node_i in range(N)]
        
        return I_total
  
class Ensemble:

  def __init__(self, No_Nodos, No_variables):
      self.No_Nodos = No_Nodos
      self.No_variables = No_variables
      self.Ensamble_original = []
      self.Ensamble_transformado = []

  def append(self, realization):
      self.Ensamble_original.append(realization)
      self.Ensamble_transformado.append(self.transform(realization))

  def average(self):
      average_realization = array(self.Ensamble_transformado[0])
      for i in range(1,len(self.Ensamble_transformado)):
          average_realization = average_realization + array(self.Ensamble_transformado[i])
      N = full_like(average_realization,1./len(self.Ensamble_transformado),float)
      return array(average_realization*N)
      
  def transform(self,realization):
      y=[]
      for i in range(self.No_Nodos):
          y.append([])
          for j in range(self.No_variables):
              y[i].append(0)

      xT=transpose(realization)
      for i in range(len(xT)):
          indice_nodo=int(i/self.No_variables)
          indice_varibale= int(i%self.No_variables)
          y[indice_nodo][indice_varibale]=array(xT[i])
      
      return y
      
class Save_Info(object):
    """ docstring for Save_Info"""
    
    def __init__(self):
        return

    def MakeDir(self,path_name):

        if not os.path.exists(path_name):
            os.makedirs(path_name)
        else:
            shutil.rmtree(path_name)
            os.makedirs(path_name)

        return

    def SaveFiles(self,path_name,make_dir=False,file_format='*.hdf5'):
    
        """
          Routine to Copy all the files with a given format into the folder path_name

          path_name: The destination folder
          format:    the files format that will be copy to path_name

        """
    
        if make_dir == True:
            self.MakeDir(path_name)
        else: pass

        #formato1 = '*.hdf5'
        #files1 = glob.iglob(os.path.join(file_format))
        #for ifile in files1: shutil.move(ifile,path_name)

        # Copy  the files in format file_format to path_name
        for ifile in glob.iglob(os.path.join(file_format)):
            # Si el archivo for formato 'file_format' ya existe, lo elimina y guarda el nuevo
            if os.path.isfile(ifile):
                try:
                    shutil.move(ifile,path_name)
                except:
                    os.remove(path_name+'/'+ifile)
                    shutil.move(ifile,path_name)

        return
    
class Measures(Models):

    def Patch_Min_HostInfected(self,X,N):
        """ RETURN THE PATCH WITH THE MAXIMUM NUMBER OF INFECTED HOST """
        A = [max(X[ith_node])  for ith_node in range(N)] 
        return A.index(min(A)) 
        
    def Patch_Min_VectorInfected(self,X,N):
        """ RETURN THE PATCH WITH THE MAXIMUM NUMBER OF INFECTED HOST """
        A = [max(X[ith_node] for ith_node in range(N))] 
        return A.index(max(A)) 
    
    def Patch_Min_TimeInf_Host(self,X,N):
        """ RETURN THE PATCH WITH THE MAXIMUM DURABILITY IN THE INFECTION """
        A = [sum(0.01*X[ith_node]  for ith_node in range(N))] 
        return A.index(max(A)) 
    
    def Patch_Min_TimeInf_Vectors(self,X,N):
        """ RETURN THE PATCH WITH THE MAXIMUM DURABILITY IN THE INFECTION """
        A = [sum(0.01*X[ith_node]  for ith_node in range(N))] 
        return A.index(max(A)) 
