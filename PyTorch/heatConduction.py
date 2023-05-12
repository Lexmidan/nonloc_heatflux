# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:15:28 2019

@author: ?
"""
import numpy as np
import pandas as pd
import utility
import time
import torch
from matplotlib import pyplot as plt
import scipy
def assemble(para, cache, alphas, betas, heatflux):
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


    ne=para['InitneProfile']
    Kb=para['boltzman']
    kappa=para['conductivity']
    dx = para['deltaX']

    numberOfNode = para['numberOfNode']
    

    #Conductivity according to notes_heatconduction.pdf
    ##Coulomb logarithm
    lamb =  np.ones(len(para['x'])) #23-np.log(np.sqrt(ne)*para['InitZbarProfile']/T**1.5)
    kappa=kappa*1.31e10/lamb*para['tau']**(betas-5/2)

    # BC informations
    typeX0 = para['x=0 type']
    valueX0 = para['x=0 value']
    typeXL = para['x=L type']
    valueXL = para['x=L value']
    
    # Containers
    T = cache['T']; T0 = cache['T0']        #let T=T[i,j] then T0=T[i, j-1]
    F = cache['F']; Jacobian = cache['Jacobian']
    
    dt = cache['dt']    
    '''    
    Loop over grid
    '''


    for i in range(0, numberOfNode):
        # BC node at x=0
        if i == 0:
            if typeX0 == 'heatFlux':
                Ug1 = utility.fixedGradient(valueX0, kappa[i], dx, T[0], alphas[i], betas[i]) #boundary values
                Jacobian[0][1] = (1/dx**2)*(alphas[i+1]*kappa[i+1]*betas[i+1]*T[i+1]**(betas[i+1]-1)*T[i]\
                                           -(betas[i+1]+1)*alphas[i+1]*kappa[i+1]*T[i+1]**betas[i+1] - alphas[i]*kappa[i]*T[i]**betas[i])
            elif typeX0 == 'fixedTemperature':
                Ug1 = utility.fixedValue(valueX0, T[1])
                Jacobian[0][1] = 0
            Jacobian[i][i] = (3/2*ne[i]*Kb)/dt + (1/dx**2)*.5*(alphas[i+1]*kappa[i+1]*T[i+1]**betas[i+1] + alphas[i]*kappa[i]*Ug1**betas[i]\
                            +2*(betas[i]+1)*alphas[i]*kappa[i]*T[i]**betas[i] - (betas[i]*alphas[i]*kappa[i]*T[i]**(betas[i]-1))*(Ug1+T[i+1]))
        # BC node at x=L
        elif i == numberOfNode-1:
            if typeXL == 'heatFlux':
                Ug2 = utility.fixedGradient(valueXL, kappa[i], dx, T[-1], alphas[i], betas[i])  #boundary values
                Jacobian[-1][-2] = (1/dx**2)*(betas[i-1]*kappa[i-1]*alphas[i-1]*T[i-1]**(betas[i-1]-1)*T[i]\
                                           -(betas[i-1]+1)*alphas[i-1]*kappa[i-1]*T[i-1]**betas[i-1] - alphas[i]*kappa[i]*T[i]**betas[i])
            elif typeXL == 'fixedTemperature':
                Ug2 = utility.fixedValue(valueXL, T[-2])
                Jacobian[-1][-2] = 0
            Jacobian[i][i] = (3/2*ne[i]*Kb)/dt + (1/dx**2)*.5*(alphas[i]*kappa[i]*Ug2**betas[i] + alphas[i-1]*kappa[i-1]*T[i-1]**betas[i-1]\
                            +2*(betas[i]+1)*alphas[i]*kappa[i]*T[i]**betas[i] - (betas[i]*alphas[i]*kappa[i]*T[i]**(betas[i]-1))*(T[i-1]+Ug2))  
        # Interior nodes

        else:   #!!! \alpha_{i+1/2} := alpha[i]
            Jacobian[i][i+1] = (1/dx**2)*.5*(alphas[i+1]*kappa[i+1]*betas[i+1]*T[i+1]**(betas[i+1]-1)*T[i]\
                                           -(betas[i+1]+1)*alphas[i+1]*kappa[i+1]*T[i+1]**betas[i+1] - alphas[i]*kappa[i]*T[i]**betas[i])
            
            Jacobian[i][i-1] = (1/dx**2)*.5*(betas[i-1]*kappa[i-1]*alphas[i-1]*T[i-1]**(betas[i-1]-1)*T[i]\
                                           -(betas[i-1]+1)*alphas[i-1]*kappa[i-1]*T[i-1]**betas[i-1] - alphas[i]*kappa[i]*T[i]**betas[i])
            
            Jacobian[i][i] = (3/2*ne[i]*Kb)/dt + (1/dx**2)*.5*(alphas[i+1]*kappa[i+1]*T[i+1]**betas[i+1] + alphas[i-1]*kappa[i-1]*T[i-1]**betas[i-1]\
                            +2*(betas[i]+1)*alphas[i]*kappa[i]*T[i]**betas[i] - (betas[i]*alphas[i]*kappa[i]*T[i]**(betas[i]-1))*(T[i-1]+T[i+1]))


    # Calculate F (right hand side vector)
    d2T = utility.secondOrder(T, Ug1, Ug2, alphas, betas,kappa)
    F = (3/2*np.array([ne]).T)*(T - T0)*Kb/dt + d2T/dx**2 # Vectorization   dT/dt - a d2T/dx2=F/dt

    # Store in cache
    cache['F'] = F; cache['Jacobian'] = Jacobian
    cache['alpha']=alphas
    cache['beta']=betas
    cache['heatflux']=heatflux
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
    T_init = para['InitTeProfile']
    T = np.reshape(T_init.values, (numberOfNode,1)) #numberOfNode rows with Tic values
    T0 = np.reshape(T_init.values, (numberOfNode,1))
    alpha_init = para['alphas']
    beta_init = para['betas']

    TProfile = np.zeros((numberOfNode, numOfTimeStep + 1))
    alpha_prof= np.zeros((numberOfNode, numOfTimeStep + 1))
    beta_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    heatflux_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    F = np.zeros((numberOfNode, 1))
    Jacobian = np.zeros((numberOfNode, numberOfNode))
    TProfile[:,0] = T.reshape(1,-1)  # first profile (column) is full of Tic 
    alpha_prof[:,0] = alpha_init.reshape(1,-1)
    beta_prof[:,0] = beta_init.reshape(1,-1)
    times=np.array([0])
    dt=Exception("dt wasn't calculated")
    kappa=Exception("kappa wasn't calculated")
    cache = {'T':T,'T0':T0,'TProfile':TProfile, 'alpha':alpha_init, 'beta':beta_init,
             'F':F,'Jacobian':Jacobian, 'time':0, 'times':times, 'dt':dt, 'conductivity': kappa,
             'Log':pd.DataFrame(), 'alpha_prof':alpha_prof, 'beta_prof':beta_prof, 'heatflux_prof':heatflux_prof}
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
    relax = para['relaxation'] 
    A = cache['Jacobian']
    B = cache['F']
    dT = np.linalg.solve(A, B)
    T = cache['T']
    T = T-dT * relax       #T(j+1)=T(j)+JI``(F)
    T[np.where(T<=0)]=10
    cache['T'] = T
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
    heatflux_prof = cache['heatflux_prof']    
    alpha = cache['alpha']      #current profile
    beta = cache['beta']
    heatflux = cache['heatflux']
    T = cache['T']
    cache['T0'] = T.copy()
    TProfile[:,timeStep] = T.reshape(1,-1)
    alpha_prof[:,timeStep] = alpha.reshape(1,-1)
    beta_prof[:,timeStep] = beta.reshape(1,-1)
    heatflux_prof[:,timeStep] = heatflux.reshape(1,-1)
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
    gradT=np.gradient(np.reshape(T, (numberOfNode)),x.values)
    #initially T had form of [[1],[2],[3]...] not [1,2,3]
    Tscaled=pd.DataFrame((np.reshape(T, (numberOfNode))-scale['T'].loc['mean'])/scale['T'].loc['std'])
    gradTscaled=pd.DataFrame((gradT-scale['gradT'].loc['mean'])/scale['gradT'].loc['std'])

    """
    Parameters given by NN:
    """    
    if para['NNmodel']==None:
        alphas=para['alphas']
        betas=para['betas']
        heatflux = para['heatflux']
    else:
        alphas, betas, heatflux = get_data_qless(para['NNmodel'], para['x'], Tscaled, gradTscaled,Zscaled, \
                                        nescaled, Knscaled, int(para['NNmodel'].fcIn.in_features/4))
                                                                                #size of the input vector
        #week 19th notes. \alpha[localtransport]:=1 ()
        #alphas[np.abs(gradT/np.reshape(T, (numberOfNode)))<1e-1]=1
        #alphas=scipy.ndimage.gaussian_filter1d(alphas,3)
        alphas[np.abs(gradT/np.reshape(T, (numberOfNode)))<1]=1
        alphas[np.abs(gradT/np.reshape(T, (numberOfNode)))>800]=20
        #alphas=scipy.ndimage.gaussian_filter1d(alphas,1)
        #alphas=scipy.ndimage.gaussian_filter1d(alphas,2)
        #betas=scipy.ndimage.gaussian_filter1d(betas,2)

    #interpolation needed in order to 'place' coefficients to the center of the cell
    ##alphas=np.interp(np.arange(0, numberOfNode)+0.5, np.arange(0,numberOfNode), alphas)
    ##betas=np.interp(np.arange(0, numberOfNode)+0.5, np.arange(0,numberOfNode), betas)
    
    cache['dt'] = para['Time_multiplier']*np.min(3/2*para['InitneProfile']*para['boltzman']*para['deltaX']**2/\
                               (para['conductivity']*para['alphas']*T[:,0]**2.5))
    cache['time']+=cache['dt']
    cache['times'] = np.append(cache['times'],cache['time'])
    log = cache['Log']
    ts = cache['ts']

    ####Iterating till resudue is small enough
    #for n in range(maxIteration):

    # cache = assemble(para, cache, alphas, betas, heatflux)
    # F = cache['F']
    # norm = np.linalg.norm(F)
    # slump = np.copy(norm)
    # n=0
    # while True:
    #     n+=1
    #     cache = assemble(para, cache, alphas, betas, heatflux)
    #     F = cache['F']
    #     norm = np.linalg.norm(F)
    #     if norm/np.mean(1.5*(np.array([para['InitneProfile']]).T)*para['boltzman']*T) < convergence:
    #         log.loc[ts,'PhysicalTime'] = cache['time']
    #         log.loc[ts,'Iteration'] = n+1
    #         log.loc[ts,'Residual'] = norm
    #         break
    #     cache = solveLinearSystem(para, cache)


    #### Limit of max number of iters
    for n in range(maxIteration):
        cache = assemble(para, cache, alphas, betas, heatflux)
        F = cache['F']
        norm = np.linalg.norm(F)
        if n==0: slump = np.copy(norm)
        if norm/np.linalg.norm(cache['T']) < convergence:
            log.loc[ts,'PhysicalTime'] = cache['time']
            log.loc[ts,'Iteration'] = n+1
            log.loc[ts,'Residual'] = norm
            break
        cache = solveLinearSystem(para, cache)
    print('[{:3.0f}'.format(ts), ']',
          '[{:6.2E}'.format(cache['time']),']',
          '[{:2.0f}'.format(n+1), ']',
          '[{:8.2E}'.format(norm/np.mean(1.5*(np.array([para['InitneProfile']]).T)*para['boltzman']*T)),']',
          '[{:8.2E}'.format(norm/slump),']',
          '[{:8.2E}'.format(np.max(cache['beta'])),']',
          '[{:8.2E}'.format(np.max(cache['alpha'])),']',
          '[{:8.2E}'.format(np.min(T)),']',
          '[{:8.2E}'.format(np.max(T)),']',
          #' [','{:8.2E}'.format(np.mean(cache['T'])),']')
          '[{:8.16E}'.format(np.mean(1.5*para['boltzman']*T*(np.array([para['InitneProfile']]).T))),']')
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
    print(' [Step] [Time] [Iter] [Residue] [Newton outcome] [Max beta] [Max alpha] [Minimal T] [Maximal T] [meanEnergy]')
    for timeStep in range(1, numOfTimeStep+1):
        cache['ts'] = timeStep
        cache = newtonIteration(para, cache)
        cache = storeUpdateResult(cache)
    TProfile = pd.DataFrame(cache['TProfile'], columns=cache['times'],index=para['x'])
    alpha_prof = pd.DataFrame(cache['alpha_prof'], columns=cache['times'],index=para['x'])
    betas_prof = pd.DataFrame(cache['beta_prof'], columns=cache['times'],index=para['x'])
    heatflux_prof = cache['heatflux_prof']
    runtime = time.time() - start
    print('[Cost] CPU time spent','%.3f'%runtime,'s')
    return TProfile, cache, alpha_prof, betas_prof, heatflux_prof


def get_data_qless(model, x, T, gradT, Z, n, Kn, lng):  
    numFields = 4 #T, gradT, Z, n, Kn
    Qdata=np.empty((0,numFields*lng), int) #2 * rad "#of points in interval" * 5 "for each phsy quantity" + 2 "for Q and beta"

    for ind, _ in enumerate(x):  #x_min=x[ind], x_max=x[ind+2*rad], x_c=x[ind+rad]
        datapoint=np.array([])          
        if ind+lng>=len(x)+1:
            break    
        else:
            datapoint=np.append(datapoint, T.iloc[ind:ind+lng]) #append all Te in xmin-xmax
            datapoint=np.append(datapoint, gradT.iloc[ind:ind+lng]) #append all gradTe in xmin-xmax
            datapoint=np.append(datapoint, Z.iloc[ind:ind+lng]) 
            datapoint=np.append(datapoint, n.iloc[ind:ind+lng]) 
            #datapoint=np.append(datapoint, Kn.iloc[ind:ind+lng]) 
            #datapoint=np.append(datapoint, x[ind:ind+lng])
            # TODO: what is the appropriate scaling here? Global (max(x)-min(x)) might be to large!
            Qdata=np.append(Qdata,[datapoint], axis=0)

    heatflux = model.heatflux_model(torch.tensor(Qdata).float()).detach().numpy()
    alphas=model.alpha_model(torch.tensor(Qdata).float()).detach().numpy()
    betas=model.beta_model(torch.tensor(Qdata).float()).detach().numpy()


    for i in range(int((len(x)-len(alphas))/2)):              
        alphas = np.append(alphas[0], alphas)
        betas = np.append(betas[0], betas)
        heatflux = np.append(heatflux[0], heatflux)

    for i in range(int((len(x)-(ind-1))/2)):              
        alphas = np.append(alphas, alphas[-1])
        betas = np.append(betas, betas[-1])
        heatflux = np.append(heatflux[0], heatflux)

    return alphas,betas,heatflux#scipy.ndimage.gaussian_filter1d(alphas,3), scipy.ndimage.gaussian_filter1d(betas,3), heatflux

#fig1, ax1 = plt.subplots(figsize=(6,3))


