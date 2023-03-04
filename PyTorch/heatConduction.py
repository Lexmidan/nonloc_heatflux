# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:15:28 2019

@author: RickFu
Modified 27.02.23
"""
import numpy as np
import pandas as pd
import parameter
import utility
import time
import torch


def assemble(para, cache):
    """ Assemble linear system Jacobian * dx = F
    
    Process:
        0. Obtain relevant informations
        1. Loop over grid:
            1.1 Deal with BC node at x=0
            1.2 Deal with BC node at x=L
            1.3 Deal with interior nodes
            1.4 Obtain values on imaginary nodes (Ug1 and Ug2)
                for two BCs
            1.4 Assemble Jacobian (a diagonal matrix)
        2. Calculate temperature gradient dT2
        3. Assemble F
    
    Return: dictionary containing cache data
    """
    
    dt = para['deltaTime']
    numberOfNode = para['numberOfNode']
    # BC informations
    typeX0 = para['x=0 type']
    valueX0 = para['x=0 value']
    typeXL = para['x=L type']
    valueXL = para['x=L value']
    
    # Containers
    T = cache['T']; T0 = cache['T0']        #let T=T[i,j] then T0=T[i, j-1]
    F = cache['F']; Jacobian = cache['Jacobian']
    
    
    x=para['x']
    #Z=para['InitZbarProfile']
    ne=para['InitneProfile']
    #Kn=para['InitKnProfile']
    scale=para['scaling']
    k=para['conductivity']
    Zscaled=para['ScaledZ']
    nescaled=para['Scaledne']
    Knscaled=para['ScaledKn']
    gradT=np.gradient(np.reshape(T, (numberOfNode)),x.values) #initially T had form of [[1],[2],[3]...] not [1,2,3]
    """
#Parameters given by NN:
    """    
        
    Tscaled=(np.reshape(T, (numberOfNode))-scale['T'].loc['mean'])/scale['T'].loc['std']
    gradTscaled=(gradT-scale['gradT'].loc['mean'])/scale['gradT'].loc['std']
    
    '''
    
    #number of points in one domain / number of features (T, gradt, Z, n, kn, x)
    lng=int(para['NNmodel'].fcIn.in_features/6)  #TODO: find a neater way to find this num

    xind=np.array([(lng-1)/2], dtype=np.uint) #store the indexes of points, where alphas and betas are calculated
    #datavector['T']=(T[:lng]-data_scaling['T'].loc['mean'])/data_scaling['T'].loc['std']
    datavector=pd.DataFrame(columns=['T','gradT', 'Z', 'n', 'Kn','x'])
    datavector['T']=Tscaled[:lng]
    datavector['gradT']=gradTscaled[:lng]
    datavector['Z']=Zscaled[:lng]
    datavector['n']=nescaled[:lng]
    datavector['Kn']=Knscaled[:lng]
    datavector['x']=np.full(lng,0) #!!! we won't need this for NN 
    alpha, beta = para['NNmodel']\
        (torch.tensor(datavector.values.flatten('F'), dtype=torch.float))
    params = pd.DataFrame([[alpha.detach().numpy(),beta.detach().numpy()]],\
                          columns=['alpha', 'beta'], index=xind)

    for ind,_ in enumerate(x, start=1):        
        if ind+lng>=len(x)+1:
            break    
        datavector=datavector.drop([0]).reset_index(drop=True)
        datavector.loc[len(datavector)+1]=[Tscaled[ind], gradTscaled[ind],\
                                           Zscaled[ind], nescaled[ind], Knscaled[ind],0]
        xind=np.append(xind, int((lng-1)/2+ind)) #index of the alphabeta
        alpha, beta = para['NNmodel']\
            (torch.tensor(datavector.values.flatten('F'), dtype=torch.float))
        params.loc[xind[-1]]=[alpha.detach().numpy(),beta.detach().numpy()]        
        if ind%5000==0:
                print(f"alpha and beta are calculated for {ind}/{len(x)-lng+1} points") 
    #Add alphas and betas at the beginning and end of intervals 
    #in order to all arrays have the same length
    #!!!!!
    for i in range(xind[0]):              
        params = pd.concat([params.iloc[0].to_frame().T.set_index(pd.Index([xind[0]-i-1])), params])
    for i in range(len(x)-xind[-1]-1):              
        params = pd.concat([params,params.iloc[0].to_frame().T.set_index(pd.Index([xind[-1]+i+1]))])
    # alphas=params['alpha']
    # betas=params['beta']
    '''
    
    alphas=np.full(len(x), 1)
    betas=np.full(len(x), 0)
    params=pd.DataFrame([alphas,betas], index=['alpha', 'beta']).T
    #!!!!!
    '''    
# Loop over grid
    '''
    for i in range(0, numberOfNode):
        
        # BC node at x=0
        if i == 0:
            dx=x[i+1]-x[i]
            temp2 = - (dt/dx**2)*((k[i]*alphas[i]+k[i+1]* alphas[i+1])\
                       *T[i]**betas[i])
            if typeX0 == 'heatFlux':
                Ug1 = utility.fixedGradient(valueX0, k, dx, T[1],i) #boundary values
                Jacobian[0][1] = temp2 * 2
            elif typeX0 == 'fixedTemperature':
                Ug1 = utility.fixedValue(valueX0, T[1])
                Jacobian[0][1] = 0
                
        # BC node at x=L
        elif i == numberOfNode-1:
            dx=x[i]-x[i-1]
            temp2 = - (dt/dx**2)*((k[i]*alphas[i]+k[i]* alphas[i])\
                       *T[i]**alphas[i])
            if typeXL == 'heatFlux':
                Ug2 = utility.fixedGradient(valueXL, k, dx, T[-2],i)  #boundary values
                Jacobian[-1][-2] = temp2 * 2
            elif typeXL == 'fixedTemperature':
                Ug2 = utility.fixedValue(valueXL, T[-2])
                Jacobian[-1][-2] = 0
                
        # Interior nodes
        else:
            dx=x[i+1]-x[i]
            temp1 =-(dt/dx**2)*(k[i] * alphas[i]*T[i-1]**betas[i-1])
            temp2 =2/(3*ne[i])+ (dt/dx**2)*(k[i]*alphas[i]+k[i+1]* alphas[i+1])\
                       *T[i]**betas[i]
            temp3 =-(dt/dx**2)*(k[i+1] * alphas[i]*T[i+1]**betas[i+1])
            Jacobian[i][i+1] = temp3
            Jacobian[i][i-1] = temp1
        Jacobian[i][i] = temp2
    
    # Calculate F (right hand side vector)
    gradq = utility.secondOrder(T, Ug1, Ug2, alphas, betas,k) #d2T/dx2
    F = (2/(3*ne[i]))*(T - T0) - (dt/dx**2)*gradq # Vectorization   dT/dt - a d2T/dx2=F/dt
    # Store in cache
    cache['F'] = -F; cache['Jacobian'] = Jacobian
    cache['alphabetas']=params
    return cache


def initialize(para):
    """ Initialize key data
    
    T: current step temperature
    T0: last step temperature
    TProfile: temperature results in time and space
    F: B as right hand side of Ax = B
    Jacobian: A as left had side of Ax = B
    
    Return: a dictionary
    """
    
    numberOfNode = para['numberOfNode']
    numOfTimeStep = para['numberOfTimeStep']
    Tic = para['InitTeProfile']
    T = np.reshape(Tic.values, (numberOfNode,1)) #numberOfNode rows with Tic values
    T0 = np.reshape(Tic.values, (numberOfNode,1))
    TProfile = np.zeros((numberOfNode, numOfTimeStep + 1))
    F = np.zeros((numberOfNode, 1))
    Jacobian = np.zeros((numberOfNode, numberOfNode))
    TProfile[:,0] = T.reshape(1,-1)  # first profile (column) is full of Tic 
    cache = {'T':T,'T0':T0,'TProfile':TProfile,
             'F':F,'Jacobian':Jacobian,
             'Log':pd.DataFrame(), 'alphabetas':pd.DataFrame()}
    return cache


def solveLinearSystem(para, cache):
    """ Solve Ax=B
    
    Process:
        1. Get A = Jacobian matrix (Jacobian)
        2. Get B = Right hand side equation (F)
        3. Calculate dT
        4. Update T
        5. Store in cache
        
    Return: a dictionary
    """
    relax = para['relaxation']   #???
    A = cache['Jacobian']
    B = cache['F']
    dT = np.linalg.solve(A, B)
    T = cache['T']
    T = dT * relax + T
    cache['T']=T
    cache['dT'] = dT
    return cache


def storeUpdateResult(cache):
    """ Store results
    Update T0
    Store temperaure results into a dataframe and 
    save it in the cache.
    """
    
    timeStep = cache['ts']
    TProfile = cache['TProfile']
    T = cache['T']
    cache['T0'] = T.copy()
    TProfile[:,timeStep] = T.reshape(1,-1)
    return cache


def newtonIteration(para, cache):
    """ Newton's Iteration for Equation System
    
    Process:
        1. Get max iteratino, convergence limit
        2. Call assemble function to get Jacobian and F(RHS)
        3. Solve for dT, update solution
        4. Evaluate F, get value of 2-norm
        5. If solution converged, break, output to screen and
           return cache.
    
    """
    
    maxIteration = para['maxIteration']
    convergence = para['convergence']
    dt = para['deltaTime']
    log = cache['Log']
    ts = cache['ts']
    for n in range(maxIteration):
        cache = assemble(para, cache)
        F = cache['F']
        norm = np.linalg.norm(F)
        if norm < convergence:
            log.loc[ts,'PhysicalTime'] = dt*ts
            log.loc[ts,'Iteration'] = n+1
            log.loc[ts,'Residual'] = norm
            break
        cache = solveLinearSystem(para, cache)
    print(' [','{:3.0f}'.format(ts), ']',
          ' [','{:6.2f}'.format(ts*dt),']',
          ' [','{:2.0f}'.format(n+1), ']',
          ' [','{:8.2E}'.format(norm),']')
    return cache


def solve(para):
    """ Main function to solve heat conduction
    
    Input: a Pandas series containing all parameters
    
    Process:
        1. Initialize cache
        2. Time marching 
        3. Newton's iteration for discretized PDE for singe time 
           step
        4. Update T, save result to T profile
    
    Return: temperature profile as final result
    """
    
    print(" Heat Conduction Solver")
    start = time.time()
    cache = initialize(para)
    numOfTimeStep = para['numberOfTimeStep']
    print(' [Step] [Pysical Time] [Iteration] [Residue]')
    for timeStep in range(1, numOfTimeStep+1):
        cache['ts'] = timeStep
        cache = newtonIteration(para, cache)
        cache = storeUpdateResult(cache)
    TProfile = cache['TProfile']
    runtime = time.time() - start
    print('[Cost] CPU time spent','%.3f'%runtime,'s')
    return TProfile, cache



if __name__ == "__main__":
    para = parameter.main()
    results, cache = solve(para)
    

