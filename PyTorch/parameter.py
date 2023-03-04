# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 22:41:22 2019

@author: RickFu
Modified on 18.02.23
"""
import postprocessing as pp
import heatConduction as hc
import pandas as pd
import numpy as np
import torch


#### Import initial profile, used in PyTorch part. (x, Te, gradTe, ne, Zbar)
init_profile=pd.read_csv('init_profile.csv', index_col=(0))


#####
init_profile=init_profile.iloc[::50,:]
init_profile.reset_index(drop=True, inplace=True)
#####

def main():
    """ Generate parameter
    
    1. Generate system-level parameters
    2. Generate material properties, grid, time, bcs
    
    Return: a Pandas series
    """
    
    column = 'values'
    df = pd.Series(name = column, dtype='float64')
    df = df.astype('object')
    
    # System-level 
    df.at['problem'] = 'HeatConduction'
    df.at['SpatialDiscretize'] = 'CenteredDifferencing'
    df.at['TimeDiscretize'] = 'BackwardEular'
    df.at['ODEsolver'] = 'NewtonIteration'
    df.at['linearSolver'] = 'numpy linalg'
    df.at['CPU'] = 1
    df.at['NNmodel']= torch.load('Model.pt')
    df.at['NNmodel'].eval()
    
    # Material
    df.at['material'] = 'steel'
    df.at['material function'] = 'constant'
    df.at['density'] = 7850
    df.at['conductivity'] = 60.5
    df.at['heatCapacity'] = 434
    
    # Grid
    df.at['length'] = 1
    df.at['numberOfNode'] = len(init_profile)
    df.at['x']=init_profile['x']
    
    # Solution
    df.at['numberOfTimeStep'] = 40#400
    df.at['deltaTime'] = 0.2
    df.at['maxIteration'] = 20
    df.at['convergence'] = 1e-2#1E-10
    df.at['relaxation'] = 1 # value in [0-1] Very sensitive!!!
    df.at['scaling']=pd.read_csv('data_scaling.csv'\
                                 , index_col=(0))
    
    # Initial conditions
    df.at['InitTeProfile'] = init_profile['Te']
    df.at['InitneProfile'] = init_profile['ne']
    df.at['InitgradTeProfile'] = init_profile['gradTe']
    df.at['InitZbarProfile'] = init_profile['Zbar']
    df.at['InitKnProfile'] = init_profile['Kn']
    
    
    # Boundary conditions
    df.at['x=0 type'] = 'heatFlux' #'heatFlux' or 'fixedTemperature'
    df.at['x=0 value'] = 0
    df.at['x=L type'] = 'heatFlux' #'heatFlux' or 'fixedTemperature'
    df.at['x=L value'] = 0  
    return df



if __name__ == "__main__":
    parameter = main()
    results, cache = hc.solve(parameter)
    T = pp.preprocess(parameter, results)
    pp.evolutionField(T)
    positions = np.linspace(parameter['x'][0], parameter['x'][-1], 10 )   #0-L  TODO: global variable?
    pp.thermalCouplePlot(T, positions)
    times = np.linspace(0, parameter['deltaTime']*parameter['numberOfTimeStep'] ,10)\
        #'numberOfTimeStep'*'deltaTime'  TODO: global variable?
    pp.temperatureDistribution(T, times)
    
    
    
    
    