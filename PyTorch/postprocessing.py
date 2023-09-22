# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 23:07:48 2019

@author: RickFu
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd


def evolutionField(results, name):
    """ Generate 3D temperature fields
    
    For better understanding of the results
    
    Inputs:
        1. parameter, a pandas series
        2. results, a numpy array
    """
    
    X = results.index
    Y = results.columns*1e9
    X, Y = np.meshgrid(X, Y)
    
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel(r'$x$ [cm]', fontsize=16,labelpad=15)
    ax.set_ylabel(r'$t$ [ns]', fontsize=16,labelpad=15)
    ax.set_zlabel(name, fontsize=16,labelpad=15)


    ax.grid(visible=None, which='minor', axis='both')
    Z = results.T.values
    ax.plot_surface(X, Y, Z, 
                    cmap=cm.seismic,
                    linewidth=0, 
                    antialiased=True)
    plt.show()



def thermalCouplePlot(results, positions):
    """ Generate x-y plots as thermo-couple data
    
    Inputs:
        1. results, a pandas DataFrame
        2. Positions, a list of positions of the generated
           grids.

    """
    
    df = results.loc[positions,:].round(2)
    df = df.T
    df = df.add_prefix('x = ')
    df = df.add_suffix(' m')
    ax = df.plot(grid=True, legend=False)
    ax.set_xlabel("Time, s")
    ax.set_ylabel("Temperature, K")



def temperatureDistribution(results, times):
    """ Generate temperature distribution at different times
    
    Inputs:
        1. results, a pandas DataFrame
        2. times, a list of timings on the calculated 
           time steps
    """
    
    df = results.loc[:,times]
    df = df.add_prefix('t = ')
    df = df.add_suffix(' s')
    ax = df.plot(grid=True)
    ax.set_xlabel("x, m")
    ax.set_ylabel("Temperature, K")



def preprocess(parameter, results, time):
    """ Pre-Process results
    
    To convert numpy array into pandas DataFrame for easier
    data processing.
    
    Input:
        1. Generated parameter serie
        2. results as a numpy array
    
    Return:
        A pandas DataFrame with index as times and 
        columns as grid positions
    """
    
    length = parameter['length']
    numberOfNode = parameter['numberOfNode']
    numOfTimeStep = parameter['numberOfTimeStep']
    grids = np.linspace(0, length, numberOfNode)#.round(5)
    times = np.linspace(0, time, numOfTimeStep+1)#.round(5)
    df = pd.DataFrame(results, 
                      index = grids, 
                      columns = times)
    return df
    
    
    
    
    
    
    
    
    
    
    
    
    