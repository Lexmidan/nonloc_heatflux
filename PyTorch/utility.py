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
            d2U[i] = (k[i] * alphas[i]/(betas[i+1]+1))*(Ug1 - 2*U[i] + U[i+1])
        elif i==(U.size - 1):
            d2U[i] = (k[i] * alphas[i]/(betas[i]+1))*(U[i-1] - 2*U[i] + Ug2)
        else:
            d2U[i] = (k[i+1] * alphas[i+1]/(betas[i+1]+1))*U[i+1]**(betas[i+1]+1)\
                -(alphas[i]*k[i-1]/(betas[i-1]+1)+alphas[i+1]*k[i+1]/(betas[i+1]+1))*U[i]**(betas[i]+1)\
                 +alphas[i-1]*k[i-1]/(betas[i-1]+1)*U[i-1]**(betas[i-1]+1)
    return d2U



