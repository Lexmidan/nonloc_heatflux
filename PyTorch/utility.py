# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:01:18 2019

@author: RickFu
Modified on 18.02.23
"""
import numpy as np


def fixedValue(value, U2):
    """  Dirichlet boundary condition

    Assume that value of variable at BC is fixed.
    Please see any numerical analysis text book for details.
    
    Return: float
    """
    
    Ug = 2 * value - U2 
    return Ug


#TEST GIT EXTS WI
def fixedGradient(q, k, dx, U1, i):
    """  Neumann boundary condition
    
    Assume that the resulted gradient at BC is fixed.
    Please see any numerical analysis text book for details.
    
    Return: float
    """
    
    Ug =  q / k[i] * 2 * dx  + U1
    return Ug



def secondOrder(U, Ug1, Ug2, alphas, betas, k):
    """ Calculate second order derivative
    
    Centered differencing approximation.
    D2U/Dx2 = (U[i-1] - 2U[i] + U[i+1])/dx**2
    
    For BC nodes, use the values on ghost nodes.
    
    Ug1: value on ghost node at x=0
    Ug2: value on ghost node at x=L
    
    Please see any numerical analysis text book for details.
    
    Return: numpy array
    """
    
    d2U = np.zeros((U.size, 1))
    for i in range(0, U.size):
            
        if i==0:
            d2U[i] = (alphas[1]*k[1]/(betas[1]+1))*U[2]**(betas[1]+1)\
                    -((alphas[0]*k[0]/(betas[0]+1))+(alphas[0]*k[1]/(betas[i]+1)))*U[1]**(betas[0]+1)\
                    +(alphas[0]*k[0]/(betas[0]+1))*U[0]**(betas[0]+1)
        elif i==(U.size - 1):
            d2U[i] = (alphas[i]*k[i]/(betas[i]+1))*Ug2**(betas[i]+1)\
                    -((alphas[i-1]*k[i-1]/(betas[i-1]+1))+(alphas[i]*k[i]/(betas[i]+1)))*U[i]**(betas[i]+1)\
                    +(alphas[i-1]*k[i-1]/(betas[i-1]+1))*U[i-1]**(betas[i-1]+1)
        else:
            d2U[i] = (alphas[i]*k[i]/(betas[i]+1))*U[i+1]**(betas[i]+1)\
                    -((alphas[i-1]*k[i-1]/(betas[i-1]+1))+(alphas[i]*k[i]/(betas[i]+1)))*U[i]**(betas[i]+1)\
                    +(alphas[i-1]*k[i-1]/(betas[i-1]+1))*U[i-1]**(betas[i-1]+1)
            
    return d2U



