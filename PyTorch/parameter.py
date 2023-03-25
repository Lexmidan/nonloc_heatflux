# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 22:41:22 2019

@author: RickFu
Modified on 18.02.23
"""
import postprocessing as pp
import heatConduction as hc
import pandas as pd
import HeatfluxModel as hfm
import numpy as np
#### Import initial profile, used in PyTorch part. (x, Te, gradTe, ne, Zbar)
init_profile=pd.read_csv('./PyTorch/init_profile.csv', index_col=(0))
#####
init_profile=init_profile.iloc[::200,:]
init_profile.reset_index(drop=True, inplace=True)
#init_profile['ne']/=init_profile['ne']*3/2
#####


def main(model):
    """ Generate parameter
    
    1. Generate system-level parameters
    2. Generate material properties, grid, time, bcs
    
    Return: a Pandas series
    """
    
    column = 'values'
    df = pd.Series(name = column, dtype='float64')
    df = df.astype('object')
    
    # System-level 
    df.at['problem'] = 'NonlocHeatConduction'
    df.at['SpatialDiscretize'] = 'CenteredDifferencing'
    df.at['TimeDiscretize'] = 'BackwardEular'
    df.at['ODEsolver'] = 'NewtonIteration'
    df.at['linearSolver'] = 'numpy linalg'
    df.at['CPU'] = 1
    
    # Grid
    df.at['length'] = init_profile['x'].iloc[-1]
    df.at['numberOfNode'] = len(init_profile)
    df.at['x']=init_profile['x']

    # Material
    df.at['material function'] = 'Given by NN'
    #df.at['conductivity']=np.ones(len(df['x']))#
    df.at['conductivity'] = (init_profile['Zbar']+0.24)/(init_profile['Zbar']*(init_profile['Zbar']+4.2))
    df.at['boltzman']=1.6e-8 #eV-> erg

    # Initial conditions
    df.at['InitTeProfile'] = init_profile['Te']
    df.at['InitneProfile'] = init_profile['ne']
    df.at['InitgradTeProfile'] = init_profile['gradTe']
    df.at['InitZbarProfile'] = init_profile['Zbar']
    df.at['InitKnProfile'] = init_profile['Kn']

    #Scaling
    df.at['scaling']=pd.read_csv('./PyTorch/data_scaling.csv', index_col=(0))
    df.at['Scaledne'] = (init_profile['ne']-df.at['scaling']['n'].loc['mean'])/df.at['scaling']['n'].loc['std']
    df.at['ScaledZ'] = (init_profile['Zbar']-df.at['scaling']['Z'].loc['mean'])/df.at['scaling']['Z'].loc['std']
    df.at['ScaledKn'] = (init_profile['Kn']-df.at['scaling']['Kn'].loc['mean'])/df.at['scaling']['Kn'].loc['std']
    
    # Boundary conditions
    df.at['x=0 type'] = 'heatFlux' #'heatFlux' or 'fixedTemperature'
    df.at['x=0 value'] = 0
    df.at['x=L type'] = 'heatFlux' #'heatFlux' or 'fixedTemperature'
    df.at['x=L value'] = 1e95
    
    #NN
    if model==None:
        df.at['NNmodel']= None
        df.at['alphas']=np.linspace(1,8, len(init_profile))#np.full(len(x), 1)
        df.at['betas']=np.full(len(init_profile), 2.5)#np.linspace(2.5,0, len(init_profile))#
    else:
        df.at['NNmodel']= model
        scale = df['scaling']
        Tscaled = (np.reshape(df['InitTeProfile'], (df['numberOfNode']))-scale['T'].loc['mean'])/scale['T'].loc['std']
        gradT = np.gradient(np.reshape(df['InitTeProfile'], (df['numberOfNode'])),df['x'].values)
        gradTscaled = (gradT-scale['gradT'].loc['mean'])/scale['gradT'].loc['std']
        df.at['alphas'], df.at['betas'] = hc.get_data_qless(df, df['x'], Tscaled, gradTscaled,df['ScaledZ'], \
                                df['Scaledne'], df['ScaledKn'], int(df['NNmodel'].fcIn.in_features/6))

    # Solution
    df.at['numberOfTimeStep'] = 70#400
    df.at['deltaX'] = df['x'].iloc[11]-df['x'].iloc[10]  #for different [i] dx differs at 16th decimal place
    df.at['deltaTime'] =np.min(3/2*df['InitneProfile']*df['boltzman']*df['deltaX']**2/\
                               (df['conductivity']*df['alphas']*df['InitTeProfile']**2.5))
    print("dt =","%8.3E" % df['deltaTime'])
    df.at['maxIteration'] = 30
    df.at['convergence'] = 1E1
    df.at['relaxation'] = 1# value in [0-1] Very sensitive!!!
    return df


#
# if __name__ == "__main__":
#     parameter = main(model)
#     results, cache, alphas, betas = hc.solve(parameter)
#     T = pp.preprocess(parameter, results)
#     pp.evolutionField(T)
#     #np.linspace(parameter['x'][0], parameter['x'].iloc[-1], 8 )   #0-L  TODO: global variable?
#     positions = T.index[::int(len(init_profile['x'])*3e-2)]
#     pp.thermalCouplePlot(T, positions)
#     times = T.columns[::int(len(T.columns)/10)][1:4]
#         #'numberOfTimeStep'*'deltaTime'  TODO: global variable?
#     pp.temperatureDistribution(T, times)
    
    
    
    