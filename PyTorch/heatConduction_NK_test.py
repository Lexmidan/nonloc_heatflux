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
import physics
import importlib
importlib.reload(scipy) 

def F_Newton_Krylov(T, para, cache, scalator):
        # BC informations

    ne = cache['ne']
    Z = cache['Zbar']
    Kn = cache['Kn']
    Kb = para['boltzman']
    x = para['x']
    dx = para['deltaX']
    dt = cache['dt']  
    T0 = cache['T0'] 

    # scaling=para['scaling']
    # Tscaled=(T-scaling['T'].loc['mean'])/scaling['T'].loc['std']
    # gradTscaled=(gradT-scaling['gradT'].loc['mean'])/scaling['gradT'].loc['std']
    # nescaled=(ne-scaling['n'].loc['mean'])/scaling['n'].loc['std']
    # Zscaled=(Z-scaling['Z'].loc['mean'])/scaling['Z'].loc['std']


    alphas, betas, heatflux = cache['alpha'], cache['beta'], cache['heatflux']
    ##Coulomb logarithm 
    coulog = np.ones(len(para['x']))  #np.ones(len(para['x'])) #23-np.log(np.sqrt(ne)*Z/T**1.5)
    ##Thermal velocity (profile)
    v=np.sqrt(T*Kb/para['m_e'])
    ##Lambda mean free path
    lamb = v**4/(ne*para['Gamma']*coulog)*1/np.sqrt(Z+1)
    gradT=np.gradient(T,x)
    ##Knudsen number accordint (5)
    Kn = -lamb*gradT/T
    kappa = para['conductivity'].values*1.31e10/coulog[:]*para['tau']**(cache['beta']-5/2)
    cache['kappa_LOCAL'] = para['conductivity']*1.31e10/coulog*para['tau']
    kQSH = 6.1e+02 # scaling constant corresponding to the SCHICK
    """
    Parameters given by NN:
    """    
    
    if para['NNmodel']==None:
        alphas = para['alphas']
        betas = para['betas']
        heatflux =  -(kQSH/Z)*((Z+0.24)/(Z+4.2))*T**2.5*gradT
    else:
        scale=para['scaling']
        alphas, betas, heatflux, cache['Kn_nonloc'] = get_data_qless(para['NNmodel'], para['x'], T, gradT, Z, \
                                        ne, Kn, int(para['NNmodel'].fcIn.in_features/4), scale)
                                                                                #size of the input vector

        # Useless for now
        dataset, ratio =qqRatio(heatflux ,cache['kappa_LOCAL'], para['x'], T, gradT, Z, \
                                                                          ne, Kn,  int(para['NNmodel'].fcIn.in_features/4))
        cache['ratio']=ratio


    F=(3/2*ne)*(T - T0)*Kb/dt + 2e6*np.gradient(heatflux, x)
    cache['alpha'], cache['beta'], cache['kappa'], cache['Kn'] = alphas, betas, kappa, Kn
    cache['heatflux']=heatflux
    cache['coulog'] = coulog
    cache['F']=F
    
    return F/scalator


def F_Newton_Krylov_SH(T, para=None, cache=None, scalator=None):
        # BC informations
    ne = cache['ne']
    Kb = para['boltzman']
    Zbar = cache['Zbar']
    x = para['x']
    dt = cache['dt']  
    T0 = cache['T0'] 
    me = 9.1094e-28 # [g]
    eV2K = 1.1604e4 # K = eV2K * eV
    erg2J = 1e-7
    kB = 1.3807e-16 * eV2K # [erg/eV]
    # Local thermal energy density
    eTh = ne * kB * T
    # Thermal velocity
    vTh = (kB * T / me)**0.5

    gradT=np.gradient(T,x)
    kQSH = 6.1e+02 # scaling constant corresponding to the SCHICK

    Qloc = -(kQSH/Zbar)*((Zbar+0.24)/(Zbar+4.2))*T**2.5*gradT
    Qfs = erg2J * vTh * eTh
    Qeff = 0.17 * Qfs * (1.0 - np.exp(-Qloc/(0.17*Qfs)))

    #WHY DA HELL IS THERE 2e6 ?! => the only way it works. np.gradient gives different result for some reason
    F=(3/2*ne)*(T - T0)*Kb/dt + 2e6*np.gradient(Qeff, x)
    return F/scalator

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

    #Define empty matrices that will contain time evolution of the profiles
    TProfile = np.zeros((numberOfNode, numOfTimeStep + 1))
    alpha_prof= np.zeros((numberOfNode, numOfTimeStep + 1))
    beta_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    heatflux_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    Zbar_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    Kn_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    ne_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    Kn_nonloc_prof= np.zeros((numberOfNode, numOfTimeStep + 1))
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
    times=np.array([0])

    coulog_init=23-np.log(np.sqrt(ne_init)*Zbar_init/T_init**1.5)

    dt=Exception("dt wasn't calculated")
    kappa=Exception("kappa wasn't calculated")
    cache = {'T':T,'T0':T0,'TProfile':TProfile, 'alpha':alpha_init, 'beta':beta_init, 'heatflux':heatflux_init,
             'F':F,'Jacobian':Jacobian, 'time':0, 'times':times, 'dt':dt, 'kappa': kappa, 'Zbar':Zbar_init, 
             'ne':ne_init,'Kn':Kn_init, 'Kn_prof':Kn_prof,'ne_prof':ne_prof,'Zbar_prof':Zbar_prof,
             'Kn_nonloc': Kn_nonloc,'Kn_nonloc_prof':Kn_nonloc_prof,'Log':pd.DataFrame(), 
             'alpha_prof':alpha_prof, 'beta_prof':beta_prof, 'heatflux_prof':heatflux_prof, 'coulog':coulog_init}
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

    alpha = cache['alpha']      #current profile
    beta = cache['beta']
    heatflux = cache['heatflux']
    Zbar = cache['Zbar']
    Kn = cache['Kn']
    ne = cache['ne']
    Kn_nonloc = cache['Kn_nonloc']
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
    return cache

def newtonIteration(para, cache, SH, verbose):
    """ Newton's Iteration for Equation System
    
    Process:
        1. Get max iteratino, convergence limit
        2. Call assemble function to get Jacobian and F(RHS)
        3. Solve for dT, update solution
        4. Evaluate F, get value of 2-norm
        5. If solution converged, break, output to screen and
           return cache.
    
    """
    

    T = cache['T'];     #let T=T[i,j] then T0=T[i, j-1] 
    #cache['dt'] = para['Time_multiplier']*np.min(3/2*para['InitneProfile']*para['boltzman']*para['deltaX']**2/((para['conductivity']*1.31e10/cache['coulog']*para['tau']**(cache['beta']-5/2))*cache['alpha']*T**2.5))
    cache['dt'] = para['Time_multiplier']*np.min(3/2*para['InitneProfile']*para['boltzman']*para['deltaX']**2/\
                               (para['conductivity']*para['alphas']*T**2.5))
    cache['time']+=cache['dt']
    cache['times'] = np.append(cache['times'],cache['time'])

    energy_init=np.mean(1.5*para['boltzman']*cache['T']*cache['ne'])
    #I use scalator because scipy takes F full of big numbers, calculates jacobian with corr-y huge numbers and assumes its inverse as zero :C
    if SH:
        firstF=F_Newton_Krylov_SH(cache['T'], para=para, cache=cache, scalator=1)
        norm = np.linalg.norm(firstF)
        F_onevar= lambda x: F_Newton_Krylov_SH(x, para=para, cache=cache, scalator=norm)  #lambda x: f(1,2,x)
    else:
        firstF=F_Newton_Krylov(cache['T'], para=para, cache=cache, scalator=1)
        norm = np.linalg.norm(firstF)
        F_onevar= lambda x: F_Newton_Krylov(x, para=para, cache=cache, scalator=norm)  #lambda x: f(1,2,x)
    

    #global numofiters
    #numofiters = 0 #counter variable for scipy


    #####IT WAS ALL ChatGPT!...mostly (Counts how many times is function called <=> number of iterations)
    def count_calls(func):
        def wrapper(*args, **kwargs):
            wrapper.calls += 1  # Increment the call count
            return func(*args, **kwargs)

        wrapper.calls = 0  # Initialize the call count
        return wrapper

    @count_calls
    def have_you_been_touched():
        return 'Yep'

    def callbackfunc(x, f, norm=norm): 
        #Function, that will return relevant info about process of jacobian solving (scipy.optimize.newton_krylov)
        #numofiters+=1
        try: 
            residue = np.abs(f/energy_init)
            slump = np.linalg.norm(f)
        except TypeError:
            slump = np.linalg.norm(f(x))
            residue = np.abs(slump/energy_init)
        have_you_been_touched()
        return residue, slump, have_you_been_touched.calls
    ####From here it was NOT ChatGPT...mostly

    #Creates an anonymous function of only one variable (as scipy requests)
    cache['T'] = scipy.optimize.newton_krylov(F_onevar, cache['T'], f_rtol=5e-9,iter=12, verbose=verbose, callback=callbackfunc)
    energy=np.mean(1.5*para['boltzman']*cache['T']*cache['ne'])
    residue, slump, numofiters= callbackfunc(cache['T'], F_onevar)
    log = cache['Log']
    ts = cache['ts']

    print('[{:3.0f}'.format(ts), ']',
          '[{:6.2E}'.format(cache['time']),']',
          '[{:2.0f}'.format(numofiters), ']',
          #'[{:8.2E}'.format(norm/energy),']',
          '[{:8.2E}'.format(residue),']', #shows what portion of energy was lost in one time step
          '[{:8.2E}'.format(slump/norm),']',
          '[{:8.2E}'.format(np.max(cache['beta'])),']',
          '[{:8.2E}'.format(np.max(cache['alpha'])),']',
          '[{:8.2E}'.format(np.min(T)),']',
          '[{:8.2E}'.format(np.max(T)),']',
          #' [','{:8.2E}'.format(np.mean(cache['T'])),']')
          '[{:8.16E}'.format(energy),']')
    cache['Log']=log #update Log
    return cache


def solve(para, SH=False,  verbose=False):
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
    print(f'Use SH ? {SH}')
    start = time.time()
    cache = initialize(para)
    numOfTimeStep = para['numberOfTimeStep']
    print(' [Step] [Time] [Iter] [Residue] [Newton outcome] [Max beta] [Max alpha] [Minimal T] [Maximal T] [meanEnergy]')
    for timeStep in range(1, numOfTimeStep+1):
        cache['ts'] = timeStep
        cache = newtonIteration(para, cache, SH, verbose)
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

            Kn_mean=np.append(Kn_mean, np.mean(Kn[ind:ind+lng])) #
            #datapoint=np.append(datapoint, Kn.iloc[ind:ind+lng]) 
            #datapoint=np.append(datapoint, x[ind:ind+lng])
            # TODO: what is the appropriate scaling here? Global (max(x)-min(x)) might be to large!
            Qdata=np.append(Qdata,[datapoint], axis=0)

    Kn_nonloc = scipy.signal.convolve(Kn, gaussian_kernel(size = lng, sigma = 6), mode='same') #length of the kernel = lng, sigma = 6, 
    heatflux = (model.forward(torch.tensor(Qdata).float())[:,0] * model.scaling['Q']['std'] + model.scaling['Q']['mean']).detach().numpy()
    beta = (model.forward(torch.tensor(Qdata).float())[:,1] * model.scaling['beta']['std'] + model.scaling['beta']['mean']).detach().numpy()
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
        Kn_mean = np.append(Kn_mean[0], Kn_mean)

    for i in range(int((len(x)-(ind-1))/2)):            
        beta = np.append(beta, beta[-1])  
        #alphas = np.append(alphas, alphas[-1])
        #betas = np.append(betas, betas[-1])
        heatflux = np.append(heatflux, heatflux[-1])
        Kn_mean = np.append(Kn_mean, Kn_mean[-1])

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



def calc_alpha(qNN, beta, Z, T, gradT, Kn):
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

    alpha = qNN/local_heatflux_beta_model
    alpha = alpha_cor(alpha, Kn)
    
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

def alpha_cor(alpha, Kn, s=2700, p=1):
    '''
    Smooth version of alpha dependence on Kn. The main idea is to grant alpha_cor=1 as Kn->0 and alpha_cor=alpha as Kn->infty
    
    args:
        alpha - 1-D (n,) np.array() The input profile alpha, independent on Kn.
        Kn - 1-D (n,) np.array() Knudsen number profile
        s, p - parameters of the smoothing function. Control the rate of alpha_cor converging to alpha as Kn->infty
    '''

    alpha_cor = 1+(s*(alpha[:]-1)*Kn[:]**(2*p))/(1+s*Kn[:]**(2*p))

    return alpha_cor

