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
from matplotlib import pyplot as plt
def assemble(para, cache, alphas, betas):
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
    #x = para['x']
    dt = para['deltaTime']
    dx = para['x'].iloc[11]-para['x'].iloc[10]     #for different [i] dx differs at 16th decimal place
    numberOfNode = para['numberOfNode']
    
    # BC informations
    typeX0 = para['x=0 type']
    valueX0 = para['x=0 value']
    typeXL = para['x=L type']
    valueXL = para['x=L value']
    
    # Containers
    T = cache['T']; T0 = cache['T0']        #let T=T[i,j] then T0=T[i, j-1]
    F = cache['F']; Jacobian = cache['Jacobian']
    
    ne=para['InitneProfile']
    Kb=1.60217e-8 #eV-> erg
    k=np.ones(len(para['x']))#para['conductivity']#
    '''    
    Loop over grid
    '''
    for i in range(0, numberOfNode):
        # BC node at x=0
        if i == 0:
            if typeX0 == 'heatFlux':
                #Ug1 = utility.fixedGradient(valueX0, k, dx, T[1],i) #boundary values
                Jacobian[0][1] = -(1/dx**2)*(k[i+1] * alphas[i+1]*T[i+1]**betas[i+1])
                F[i] = (3/2*ne[i])*(T[i] - T0[i])*Kb/dt - (1/dx**2)*(k[i]*alphas[i])/(betas[i]+1)*(T[i+1]**(betas[i+1] + 1.0)-T[i]**(betas[i] + 1.0))
            elif typeX0 == 'fixedTemperature':
                #Ug1 = utility.fixedValue(valueX0, T[1])
                F[i]=0
                Jacobian[0][1] = 0
            Jacobian[i][i] = (3/2*ne[i]*Kb)/dt+ (1/dx**2)*(2*k[i]* alphas[i])\
                            *T[i]**betas[i]
        # BC node at x=L
        elif i == numberOfNode-1:
            if typeXL == 'heatFlux':
                #Ug2 = utility.fixedGradient(valueXL, k, dx, T[-2],i)  #boundary values
                F[i]=(3/2*ne[i])*(T[i] - T0[i])*Kb/dt - (1/dx**2)*(k[i]*alphas[i])/(betas[i]+1)*(T[i-1]**(betas[i-1] + 1.0)-T[i]**(betas[i] + 1.0))
                Jacobian[-1][-2] = -(1/dx**2)*(k[i-1] * alphas[i-1]*T[i-1]**betas[i-1])
            elif typeXL == 'fixedTemperature':
                #Ug2 = utility.fixedValue(valueXL, T[-2])
                F[i]=0
                Jacobian[-1][-2] = 0
            Jacobian[i][i] = (3/2*ne[i]*Kb)/dt + (1/dx**2)*(k[i-1]*alphas[i-1] + k[i]* alphas[i])\
                            *T[i]**betas[i]  
        # Interior nodes
        else:   #!!! \alpha_{i+1/2} := alpha[i]
            Jacobian[i][i+1] = -(1/dx**2)*(k[i+1] * alphas[i+1]*T[i+1]**betas[i+1])
            Jacobian[i][i-1] = -(1/dx**2)*(k[i-1] * alphas[i-1]*T[i-1]**betas[i-1])
            Jacobian[i][i] = (3/2*ne[i]*Kb)/dt + (1/dx**2)*(k[i-1]*alphas[i-1] + k[i]* alphas[i])\
                           *T[i]**betas[i]
            F[i] = (3/2*ne[i])*(T[i] - T0[i])*Kb/dt - ((alphas[i]*k[i]/(betas[i]+1))*T[i+1]**(betas[i]+1)\
                    -((alphas[i-1]*k[i-1]/(betas[i-1]+1)) + (alphas[i]*k[i]/(betas[i]+1)))*T[i]**(betas[i]+1)\
                    +(alphas[i-1]*k[i-1]/(betas[i-1]+1))*T[i-1]**(betas[i-1]+1))/dx**2 # Vectorization   dT/dt - a d2T/dx2=F/dt
    
    # Calculate F (right hand side vector)
    #d2T = utility.secondOrder(T, Ug1, Ug2, alphas, betas,k) #d2T/dx2
    #F = (3/2*np.array([ne]).T)*(T - T0)*Kb/dt - d2T/dx**2 # Vectorization   dT/dt - a d2T/dx2=F/dt
    # Store in cache
    cache['F'] = F; cache['Jacobian'] = Jacobian
    cache['alpha']=alphas
    cache['beta']=betas
    #print('it num')
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
    alpha = np.ones((numberOfNode))
    beta = np.zeros((numberOfNode))
    TProfile = np.zeros((numberOfNode, numOfTimeStep + 1))
    alpha_prof= np.zeros((numberOfNode, numOfTimeStep + 1))
    beta_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    F = np.zeros((numberOfNode, 1))
    Jacobian = np.zeros((numberOfNode, numberOfNode))
    TProfile[:,0] = T.reshape(1,-1)  # first profile (column) is full of Tic 
    cache = {'T':T,'T0':T0,'TProfile':TProfile, 'alpha':alpha, 'beta':beta,
             'F':F,'Jacobian':Jacobian,
             'Log':pd.DataFrame(), 'alpha_prof':alpha_prof, 'beta_prof':beta_prof,  }
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
    T = T-dT * relax       #T(j+1)=T(j)+JI``(F)
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
    alpha_prof = cache['alpha_prof']     #all profiles of params
    beta_prof = cache['beta_prof']       
    alpha = cache['alpha']      #current profile
    beta = cache['beta']
    T = cache['T']
    cache['T0'] = T.copy()
    TProfile[:,timeStep] = T.reshape(1,-1)
    alpha_prof[:,timeStep] = alpha.reshape(1,-1)
    beta_prof[:,timeStep] = beta.reshape(1,-1)
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
    
    '''
    NN PART
    '''
    
    x = para['x']
    numberOfNode = para['numberOfNode']
    T = cache['T'];     #let T=T[i,j] then T0=T[i, j-1]
    scale=para['scaling']
    Zscaled=para['ScaledZ']
    nescaled=para['Scaledne']
    Knscaled=para['ScaledKn']
    gradT=np.gradient(np.reshape(T, (numberOfNode)),x.values) #initially T had form of [[1],[2],[3]...] not [1,2,3]
    Tscaled=(np.reshape(T, (numberOfNode))-scale['T'].loc['mean'])/scale['T'].loc['std']
    gradTscaled=(gradT-scale['gradT'].loc['mean'])/scale['gradT'].loc['std']
    """
    Parameters given by NN:
    """    
    # alphas, betas = get_data_qless(para, para['x'], Tscaled, gradTscaled,Zscaled, \
    #                                 nescaled, Knscaled, int(para['NNmodel'].fcIn.in_features/6)) #number of points
    # params=pd.DataFrame([alphas,betas], index=['alpha', 'beta']).T
    
    
    alphas=np.linspace(1,8, len(x))#np.full(len(x), 1)
    betas=np.linspace(2.5,0, len(x))#np.full(len(x), 2.5)
    alphas=np.interp(np.arange(0, numberOfNode)+0.5, np.arange(0,numberOfNode), alphas)
    betas=np.interp(np.arange(0, numberOfNode)+0.5, np.arange(0,numberOfNode), betas)
    
    
    dt = para['deltaTime']
    log = cache['Log']
    ts = cache['ts']
    for n in range(maxIteration):
        cache = assemble(para, cache, alphas, betas)
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
          ' [','{:8.2E}'.format(norm),']',
          ' [','{:8.2E}'.format(np.min(cache['T'])),']',
          ' [','{:8.2E}'.format(np.max(cache['T'])),']',
          #' [','{:8.2E}'.format(np.mean(cache['T'])),']')
          ' [','{:8.2E}'.format(np.mean(cache['T']*(np.array([para['InitneProfile']]).T))),']')
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
    print(' [Step] [Pysical Time] [Iteration] [Residue] [minT] [maxT] [E]')
    for timeStep in range(1, numOfTimeStep+1):
        cache['ts'] = timeStep
        cache = newtonIteration(para, cache)
        
        cache = storeUpdateResult(cache)
    TProfile = cache['TProfile']
    alpha_prof = cache['alpha_prof']
    betas_prof = cache['beta_prof']
    runtime = time.time() - start
    print('[Cost] CPU time spent','%.3f'%runtime,'s')
    return TProfile, cache, alpha_prof, betas_prof

def get_data_qless(para, x, T, gradT, Z, n, Kn, lng):  
    numFields = 5 #T, gradT, Z, n, Kn
    Qdata=np.empty((0,numFields*lng), int) #2 * rad "#of points in interval" * 5 "for each phsy quantity" + 2 "for Q and beta"
    for ind, _ in enumerate(x):  #x_min=x[ind], x_max=x[ind+2*rad], x_c=x[ind+rad]
        datapoint=np.array([])          
        if ind+lng>=len(x)+1:
            break    
        else:
            datapoint=np.append(datapoint, T[ind:ind+lng]) #append all Te in xmin-xmax
            datapoint=np.append(datapoint, gradT[ind:ind+lng]) #append all gradTe in xmin-xmax
            datapoint=np.append(datapoint, Z[ind:ind+lng]) #append all Zbar in xmin-xmax
            datapoint=np.append(datapoint, n[ind:ind+lng]) #append all gradTe in xmin-xmax
            datapoint=np.append(datapoint, Kn[ind:ind+lng]) #append all Knudsen number in xmin-xmax
            # TODO: what is the appropriate scaling here? Global (max(x)-min(x)) might be to large!
            Qdata=np.append(Qdata,[datapoint], axis=0)
            # if ind%5000==0:
            #     print(f"We're done with {ind}/{len(x)-lng+1} points") 
    Qdata =torch.tensor(np.c_[Qdata, np.empty([len(Qdata), lng])]).float()
    
    heatflux = para['NNmodel'].heatflux_model(Qdata.float()).detach().numpy()
    ax1.plot(heatflux)
    alphas=para['NNmodel'].alpha_model(Qdata.float()).detach().numpy()
    betas=para['NNmodel'].beta_model(Qdata.float()).detach().numpy()
    for i in range(int((len(x)-len(alphas))/2)):              
        alphas = np.append(alphas[0], alphas)
        betas = np.append(betas[0], betas)
    for i in range(int((len(x)-(ind-1))/2)):              
        alphas = np.append(alphas, alphas[-1])
        betas = np.append(betas, betas[-1])
    return alphas, betas

fig1, ax1 = plt.subplots(figsize=(6,3))


if __name__ == "__main__":
    para = parameter.main()
    results, cache, alphas, betas = solve(para)
    

