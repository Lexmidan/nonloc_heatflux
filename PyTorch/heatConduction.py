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


    # BC informations
    typeX0 = para['x=0 type']
    valueX0 = para['x=0 value']
    typeXL = para['x=L type']
    valueXL = para['x=L value']
    ne = cache['ne']
    Z = cache['Zbar']
    Kn = cache['Kn']
    Kb = para['boltzman']
    x = para['x']
    dx = para['deltaX']

    numberOfNode = para['numberOfNode']
    # Containers
    T = cache['T']; T0 = cache['T0']        #let T=T[i,j] then T0=T[i, j-1]
    F = cache['F']; Jacobian = cache['Jacobian']
    dt = cache['dt']  


    alphas, betas, heatflux = cache['alpha'], cache['beta'], cache['heatflux']
        ##Coulomb logarithm 
    coulog = 23-np.log(np.sqrt(ne)*Z/T**1.5) #np.ones(len(para['x'])) #23-np.log(np.sqrt(ne)*Z/T**1.5)
    
        ##Thermal velocity (profile)
    v=np.sqrt(T*Kb/para['m_e'])
        ##Lambda mean free path
    lamb = v**4/(ne*para['Gamma']*coulog)*1/np.sqrt(Z+1)
    gradT=np.gradient(T,x)
    ##Knudsen number accordint (5)
    Kn = -lamb*gradT/T
    kappa = para['conductivity']*1.31e10/coulog*para['tau']**(cache['beta']-5/2)
    cache['kappa_LOCAL'] = para['conductivity']*1.31e10/coulog*para['tau']
    kQSH = 6.1e+02
    #Kn= np.sqrt(T*Kb/para['m_e'])**4/(ne*(4 * const.pi * df['q_e']**4/df['m_e']**2)*(23-np.log(np.sqrt(ne)*Z/T**1.5)))*1/np.sqrt(Z+1)*gradT/T

    """
    Parameters given by NN:
    """    
    
    if para['NNmodel']==None:
        alphas = para['alphas']
        betas = para['betas']
        heatflux = -(kQSH/Z)*((Z+0.24)/(Z+4.2))*T**2.5*gradT
        if cache['FluxLimiter']:
            # Local thermal energy density
            erg2J = 1e-7
            eTh = ne * Kb * T
            Qfs = erg2J * v * eTh
            heatflux = 0.17 * Qfs * (1.0 - np.exp(-heatflux/(0.17*Qfs)))
            Kn_nonloc = scipy.signal.convolve(Kn, gaussian_kernel(size = 23, sigma = 6), mode='same')
            alphas = calc_alpha(heatflux, betas, Z, T, gradT, Kn_nonloc, AdjustAlpha=False)
    elif para['NNeachiter'] and not cache['FluxLimiter']:
        scale=para['scaling']
        alphas, betas, heatflux, cache['Kn_nonloc'] = get_data_qless(para['NNmodel'], para['x'], T, gradT, Z, \
                                        ne, Kn, int(para['NNmodel'].fcIn.in_features/4), scale)
                                                                                #size of the input vector

        dataset, qqratio =qqRatio(heatflux ,cache['kappa_LOCAL'], para['x'], T, gradT, Z, \
                                                                          ne, Kn,  int(para['NNmodel'].fcIn.in_features/4))
        cache['qqratio']=qqratio.values
    elif para['NNmodel']!=None and cache['FluxLimiter']:
        raise Exception('Cannot use both model and limiter')


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
    F = (3/2*ne)*(T - T0)*Kb/dt + d2T/dx**2 # Vectorization   dT/dt - a d2T/dx2=F/dt

    # Store in cache
    cache['F'] = F; cache['Jacobian'] = Jacobian
    cache['alpha'], cache['beta'], cache['kappa'], cache['Kn'], cache['heatflux'] = alphas, betas, kappa, Kn, heatflux
    cache['coulog'] = coulog
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
    Zbar_init=para['InitZbarProfile']
    ne_init=para['InitneProfile']
    Kn_init=para['InitKnProfile']
    Kn_nonloc=para['Kn_nonloc']
    T = T_init 
    T0 = T_init
    alpha_init = para['alphas']
    beta_init = para['betas']
    heatflux_init= para['heatflux']
    qqratio_init = np.zeros(len(T))

    #Define empty matrices that will contain time evolution of the profiles
    TProfile = np.zeros((numberOfNode, numOfTimeStep + 1))
    alpha_prof= np.zeros((numberOfNode, numOfTimeStep + 1))
    beta_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    heatflux_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    Zbar_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    Kn_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    ne_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    Kn_nonloc_prof= np.zeros((numberOfNode, numOfTimeStep + 1))
    qqratio_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    F = np.zeros((numberOfNode, 1))
    Jacobian = np.zeros((numberOfNode, numberOfNode))

    #Filling first column with initial values of the quantities
    TProfile[:,0] = T.reshape(1,-1)
    alpha_prof[:,0] = alpha_init.reshape(1,-1)
    beta_prof[:,0] = beta_init.reshape(1,-1)
    heatflux_prof[:,0] = heatflux_init.reshape(1,-1)
    Kn_prof[:,0] = Kn_init.reshape(1,-1)
    ne_prof[:,0] = ne_init.reshape(1,-1)
    Zbar_prof[:,0] = Zbar_init.reshape(1,-1)
    qqratio_prof[:,0] = qqratio_init.reshape(1,-1)
    times=np.array([0])

    coulog_init=23-np.log(np.sqrt(ne_init)*Zbar_init/T_init**1.5)

    dt=Exception("dt wasn't calculated")
    kappa=Exception("kappa wasn't calculated")
    cache = {'T':T,'T0':T0,'TProfile':TProfile, 'alpha':alpha_init, 'beta':beta_init, 'heatflux':heatflux_init,
             'F':F,'Jacobian':Jacobian, 'time':0, 'times':times, 'dt':dt, 'kappa': kappa, 'Zbar':Zbar_init, 
             'ne':ne_init,'Kn':Kn_init, 'Kn_prof':Kn_prof,'ne_prof':ne_prof,'Zbar_prof':Zbar_prof,
             'Kn_nonloc': Kn_nonloc,'Kn_nonloc_prof':Kn_nonloc_prof,'Log':pd.DataFrame(), 'qqratio': qqratio_init, 'qqratio_prof':qqratio_prof,
             'alpha_prof':alpha_prof, 'beta_prof':beta_prof, 'heatflux_prof':heatflux_prof, 'coulog':coulog_init}
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
    Zbar_prof = cache['Zbar_prof']   
    Kn_prof = cache['Kn_prof']   
    ne_prof = cache['ne_prof']
    Kn_nonloc_prof = cache['Kn_nonloc_prof']
    qqratio_prof = cache['qqratio_prof']


    alpha = cache['alpha']      #current profile
    beta = cache['beta']
    heatflux = cache['heatflux']
    Zbar = cache['Zbar']
    Kn = cache['Kn']
    ne = cache['ne']
    Kn_nonloc = cache['Kn_nonloc']
    qqratio=cache['qqratio']
    T = cache['T']
    cache['T0'] = T.copy()

    TProfile[:,timeStep] = T.reshape(1,-1)
    alpha_prof[:,timeStep] = alpha.reshape(1,-1)
    beta_prof[:,timeStep] = beta.reshape(1,-1)
    heatflux_prof[:,timeStep] = heatflux.reshape(1,-1)
    Zbar_prof[:,timeStep] = Zbar.reshape(1,-1)
    Kn_prof[:,timeStep] = Kn.reshape(1,-1)
    ne_prof[:,timeStep] = ne.reshape(1,-1)
    Kn_nonloc_prof[:,timeStep] = Kn_nonloc.reshape(1,-1)
    qqratio_prof[:,timeStep] = qqratio.reshape(1,-1)
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

    T = cache['T'];     #let T=T[i,j] then T0=T[i, j-1] 
    #cache['dt'] = para['Time_multiplier']*np.min(3/2*para['InitneProfile']*para['boltzman']*para['deltaX']**2/((para['conductivity']*1.31e10/cache['coulog']*para['tau']**(cache['beta']-5/2))*cache['alpha']*T**2.5))
    cache['dt'] = para['Time_multiplier']*np.min(3/2*para['InitneProfile']*para['boltzman']*para['deltaX']**2/\
                               (para['conductivity']*para['alphas']*T**2.5))
    cache['time']+=cache['dt']
    cache['times'] = np.append(cache['times'],cache['time'])
    log = cache['Log']
    ts = cache['ts']


    if ~para['NNeachiter'] and para['NNmodel']!=None:
        ne = cache['ne']
        Z = cache['Zbar']
        Kn = cache['Kn']
        x = para['x']
        gradT=np.gradient(T,x)
        coulog = 23-np.log(np.sqrt(ne)*Z/T**1.5) #np.ones(len(para['x'])) #23-np.log(np.sqrt(ne)*Z/T**1.5)
        Kb = para['boltzman']
            ##Thermal velocity (profile)
        v=np.sqrt(T*Kb/para['m_e'])
            ##Lambda mean free path
        lamb = v**4/(ne*para['Gamma']*coulog)*1/np.sqrt(Z+1)
        gradT=np.gradient(T,x)
        ##Knudsen number accordint (5)
        Kn = -lamb*gradT/T
        scale=para['scaling']
        cache['alpha'], cache['beta'], cache['heatflux'], cache['Kn_nonloc'] = get_data_qless(para['NNmodel'], para['x'], T, gradT, Z, \
                                        ne, Kn, int(para['NNmodel'].fcIn.in_features/4), scale)
                                                                                #size of the input vector




    if para['Break_condition']=='max_iter':
        for n in range(maxIteration):
            cache = assemble(para, cache)
            F = cache['F']
            norm = np.linalg.norm(F)
            energy=np.mean(1.5*para['boltzman']*cache['T']*cache['ne'])
            
            if n==0: slump, energy_init = np.copy(norm), np.copy(energy)
            #if norm/np.linalg.norm(cache['T']) < convergence:
            elif np.abs(energy_init-energy)/energy_init < convergence and n!=0:
                log.loc[ts,'PhysicalTime'] = cache['time']
                log.loc[ts,'Iteration'] = n+1
                log.loc[ts,'Residual'] = np.abs(energy_init-energy)/energy_init
                break
            cache = solveLinearSystem(para, cache)

    elif para['Break_condition']=='lower_bound':
        n=0
        while True:
            cache = assemble(para, cache)
            F = cache['F']
            norm = np.linalg.norm(F)
            energy=np.mean(1.5*para['boltzman']*cache['T']*cache['ne'])
            if n==0: slump, energy_init = np.copy(norm), np.copy(energy)
            cache = solveLinearSystem(para, cache)
            n+=1

            if np.abs(energy_init-energy)/energy_init < convergence:
                break
    else: 
        print('Wrong break condition')
        quit()

    print('[{:3.0f}'.format(ts), ']',
          '[{:6.2E}'.format(cache['time']),']',
          '[{:2.0f}'.format(n+1), ']',
          #'[{:8.2E}'.format(norm/energy),']',
          '[{:8.2E}'.format(np.abs(energy_init-energy)/energy_init),']', #shows what portion of energy was lost in one time step
          '[{:8.2E}'.format(norm/slump),']',
          '[{:8.2E}'.format(np.max(cache['beta'])),']',
          '[{:8.2E}'.format(np.max(cache['alpha'])),']',
          '[{:8.2E}'.format(np.min(T)),']',
          '[{:8.2E}'.format(np.max(T)),']',
          #' [','{:8.2E}'.format(np.mean(cache['T'])),']')
          '[{:8.16E}'.format(energy),']')
    cache['Log']=log #update Log
    return cache


def solve(para, FluxLimiter=False):
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
    cache['FluxLimiter']=FluxLimiter
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


def get_data_qless(model, x, T, gradT, Z, n, Kn, lng, scaling):  
    
    """
    Takes model and profiles of physical values, scales them and compiles data frame, where every row corresponds to input vector of NN,
    and every next row is assembled of interval moved a one point to the right relativily to the previous row (sliding interval).
    Calculates non local heatflux with NN model, returns that heatflux. 

    args: 

        model - Trained and evaluated PyTorch model
        x - np.array() cont. coordinates of every point of the profile, basically is used only for indexing
        T, gradT, Z, n - np.arrays() of SCALED values of the physical quantities, which ARE used in NN
        KnUnscaled - np.array() of UNSCALED Knudsen number, needed only for comparing local and non local heatfluxes
        lng - number of points used for every quantity. Defines how 

    output:
        heatflux - np.array() of heatfluxes calculated with NN. Number of elements is equal to len(x), 
                    BUT roughly lng/2 of first and last points are equal due to the last two 'for' cycles
    """
                #alphas betas - calculated fully from NN. alpha and beta are calculated manually

    Tscaled=(T-scaling['T'].loc['mean'])/scaling['T'].loc['std']
    gradTscaled=(gradT-scaling['gradT'].loc['mean'])/scaling['gradT'].loc['std']
    nescaled=(n-scaling['n'].loc['mean'])/scaling['n'].loc['std']
    Zscaled=(Z-scaling['Z'].loc['mean'])/scaling['Z'].loc['std']

    Kn_mean=np.array([])  #Will contain mean values of Kn for each interval
    nonloc_tester=pd.DataFrame(columns=['xmin', 'xmax', 'is loc']) #Will contain information whether the Kn<10e-3 is complied for the whole vector
    numFields = 4 #T, gradT, Z, n
    Qdata=np.empty((0,numFields*lng), int) #2 * rad "#of points in interval" * 5 "for each phsy quantity" + 2 "for Q and beta" 
    
    for ind, _ in enumerate(x):  #x_min=x[ind], x_max=x[ind+2*rad], x_c=x[ind+rad]
        datapoint=np.array([])  

        if ind+lng>=len(x):
            break    
        else:
            datapoint=np.append(datapoint, Tscaled[ind:ind+lng]) #append all Te in xmin-xmax
            datapoint=np.append(datapoint, gradTscaled[ind:ind+lng]) #append all gradTe in xmin-xmax
            datapoint=np.append(datapoint, Zscaled[ind:ind+lng]) 
            datapoint=np.append(datapoint, nescaled[ind:ind+lng]) 
            nonloc_tester.loc[len(nonloc_tester.index)]=[x[ind],x[ind+lng], any(np.abs(Kn[ind:ind+lng])<1e-3)]

            #Kn_mean=np.append(Kn_mean, np.mean(Kn[ind:ind+lng])) #
            #datapoint=np.append(datapoint, Kn.iloc[ind:ind+lng]) 
            #datapoint=np.append(datapoint, x[ind:ind+lng])
            # TODO: what is the appropriate scaling here? Global (max(x)-min(x)) might be to large!
            Qdata=np.append(Qdata,[datapoint], axis=0)

    Kn_nonloc = scipy.signal.convolve(Kn, gaussian_kernel(size = lng, sigma = 6), mode='same') #length of the kernel = lng, sigma = 6, 
    heatflux = (model.forward(torch.tensor(Qdata).float())[:,0] * model.scaling['Q']['std'] + model.scaling['Q']['mean']).detach().numpy()
    #heatflux = scipy.ndimage.gaussian_filter(heatflux, sigma=2)
    #heatflux = scipy.signal.savgol_filter(heatflux, 11, 2)
    beta = (model.forward(torch.tensor(Qdata).float())[:,1] * model.scaling['beta']['std'] + model.scaling['beta']['mean']).detach().numpy()
    #beta = scipy.ndimage.gaussian_filter(beta, sigma=2)
    #beta = scipy.signal.savgol_filter(beta, 11, 2)
    beta[beta<1e-6] = 1e-6 # Make sure the power of diffusivity is positive

    #####OLD INTERFACE using AlphaBetaModel

    #alphas=model.alpha_model(torch.tensor(Qdata).float()).detach().numpy() 
    #betas=model.beta_model(torch.tensor(Qdata).float()).detach().numpy()

    ######
    

    #appending values to the end and to the beginning of the profile in order to alpha, 
    #beta, heatflux have all the same size corresponding to x
    for i in range(int((len(x)-len(beta))/2)):              
        beta = np.append(beta[0], beta)
        #alphas = np.append(alphas[0], alphas)
        #betas = np.append(betas[0], betas)
        heatflux = np.append(heatflux[0], heatflux)
        #Kn_mean = np.append(Kn_mean[0], Kn_mean)

    for i in range(int((len(x)-(ind-1))/2)):            
        beta = np.append(beta, beta[-1])  
        #alphas = np.append(alphas, alphas[-1])
        #betas = np.append(betas, betas[-1])
        heatflux = np.append(heatflux, heatflux[-1])
        #Kn_mean = np.append(Kn_mean, Kn_mean[-1])
    #Flux limiter

    #alpha calculated with 
    alpha = calc_alpha(heatflux, beta, Z, T, gradT, Kn_nonloc) #Kn_mean

    return alpha, beta, heatflux, Kn_nonloc #Kn_mean




def qqRatio(qNN, kappa, x, T, gradT, Z, n, KnUnscaled, lng):  
    """
    Takes nonlocal heatflux calculated with NN and values of physical quantities needed for calculaing the local heatflux, 
    compares them, returns ratio and DataFrame full of ratios for every *sliding* interval of length=lng, in corresp. with get_data_qless()

    args:
        qNN - (n,) np.array(), heatflux calculated with NN
        kappa, T, gradT - (n,) np.array(), from which local heatflux will be calculated
        x - (n,) np.array(), used basically only for indexing the sliding interval
    """
    qloc=kappa * T**2.5 * gradT
    Ratio = qNN/qloc
    Qdata=np.empty((0,lng+1), int) #2 * rad "#of points in interval" * 5 "for each phsy quantity" + 2 "for Q and beta"

    for ind, _ in enumerate(x):  #x_min=x[ind], x_max=x[ind+2*rad], x_c=x[ind+rad]
        datapoint=np.array([])          
        if ind+lng>=len(x):
            break    
        else:
            datapoint=np.append(datapoint, ind)
            datapoint=np.append(datapoint, Ratio[ind:ind+lng])
            Qdata=np.append(Qdata,[datapoint], axis=0)

    df=pd.DataFrame(Qdata[:,1:], index=Qdata[:,0]) #Data is intervals of ratio, index corresponds to where the interval starts on x axis.

    return df, Ratio


def calc_alpha(qNN, beta, Z, T, gradT, Kn, AdjustAlpha=True):
    '''
    Calculates local heatflux, then by comparing the latter with heatflux given by NN calculates alpha, after which 
    alpha is adjusted to the principle (Kn~0 => alpha:=1)

    args:
        *all - (n,) np.array()

    output:
        alpha - (n,) np.array()
    '''

    # Get beta local heat flux
    TkeV = 1e-3 * T[:]; gradTkeV = 1e-3 * gradT[:]
    kQSH = 6.1e+02 * 1e3**2.5 * 1e3 # scaling constant consistent with SCHICK and T in keV
    local_heatflux_beta_model = - kQSH / Z[:] * ((Z[:] + 0.24) / (Z[:] + 4.2))\
      * TkeV[:]**beta * gradTkeV[:]
    
    local_heatflux_beta_model[local_heatflux_beta_model<1e-3]=qNN[local_heatflux_beta_model<1e-3]

    alpha = qNN/local_heatflux_beta_model

    if AdjustAlpha:
        alpha = alpha_cor(alpha, Kn)
        #alpha = alpha_cor2(alpha, np.arange(len(alpha)))#Kn)
    
    #alpha[alpha>100]=1
    return alpha  
#fig1, ax1 = plt.subplots(figsize=(6,3))


def gaussian_kernel(size,sigma):
    '''
    Gauss kernel of given size

    args: 
        size - int.
        sigma - int. Standard deviation in gaussian.
    '''
    filter_range = np.linspace(-int(size/2),int(size/2),size)
    gaussian_filter = [1 / (sigma * np.sqrt(2*np.pi)) * np.exp(-x**2/(2*sigma**2)) for x in filter_range]
    return gaussian_filter

def alpha_cor(alpha, Kn, s=2.5e5, p=1):
    '''
    Smooth version of alpha dependence on Kn. The main idea is to grant alpha_cor=1 as Kn->0 and alpha_cor=alpha as Kn->infty
    
    args:
        alpha - 1-D (n,) np.array() The input profile alpha, independent on Kn.
        Kn - 1-D (n,) np.array() Knudsen number profile
        s, p - parameters of the smoothing function. Control the rate of alpha_cor converging to alpha as Kn->infty
    '''

    alpha_cor = 1+(s*(alpha[:]-1)*Kn[:]**(2*p))/(1+s*Kn[:]**(2*p))

    return alpha_cor

def alpha_cor2(alpha, x, s=1e-2, p=130, d=380):
    '''
    Smooth version of alpha dependence on Kn. The main idea is to grant alpha_cor=1 as Kn->0 and alpha_cor=alpha as Kn->infty
    
    args:
        alpha - 1-D (n,) np.array() The input profile alpha, independent on Kn.
        Kn - 1-D (n,) np.array() Knudsen number profile
        s, p - parameters of the smoothing function. Control the rate of alpha_cor converging to alpha as Kn->infty
    '''
    alpha[:p]=1
    alpha_cor = 1+(s*(alpha[:]-1)*(x[:]-p)**(2))/(1+s*(x[:]-p)**(2))
    alpha_cor = 1+(s*(alpha[:]-1)*(x[:]-d)**(2))/(1+s*(x[:]-d)**(2))
    alpha_cor[d:]=1
    return alpha_cor

def gaussian_kernel(size,sigma):
    '''
    Gauss kernel of given size

    args: 
        size - int.
        sigma - int. Standard deviation in gaussian.
    '''
    filter_range = np.linspace(-int(size/2),int(size/2),size)
    gaussian_filter = [1 / (sigma * np.sqrt(2*np.pi)) * np.exp(-x**2/(2*sigma**2)) for x in filter_range]
    return gaussian_filter

