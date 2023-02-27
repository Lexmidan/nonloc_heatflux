# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 10:53:26 2023

@author: aleks
"""



import numpy as np
import pandas as pd
import torch

init_profile=pd.read_csv('C:/Users/aleks/Projects/ML-student-projects/PyTorch/init_profile.csv', index_col=(0))
data_scaling=pd.read_csv('C:/Users/aleks/Projects/ML-student-projects/PyTorch/data_scaling.csv', index_col=(0))

    #### If we're OK with much fewer alphabetas, distributed equidistantly along xaxis 
    # alphas=np.array([])
    # betas=np.array([])
    # xind=np.array([], dtype=np.uint) #store the indexes of points, where alphas and betas are calculated
    # j=0
    # while int(j+rad)<=len(x):
    #     datavector=np.array([])
    #     datavector=np.append(datavector, T[j:j+rad])
    #     datavector=np.append(datavector, gradT[j:j+rad])
    #     datavector=np.append(datavector, Z[j:j+rad])
    #     datavector=np.append(datavector, ne[j:j+rad])
    #     datavector=np.append(datavector, np.full(21,0))  # Knx we won't need this for NN 
    #     datavector=np.append(datavector, np.full(21,0))  # x
    #     alpha, beta = para['NNmodel'](torch.tensor(datavector, dtype=torch.float))
    #     xind=np.append(xind, 11+j) #index of the center of alpha/beta pair in xref frame
    #     alphas=np.append(alphas,alpha.detach().numpy())
    #     betas=np.append(betas,beta.detach().numpy())
    #     j+=rad
    #     params = list(zip(alphas, betas))
    # params=pd.DataFrame(params,columns=['alpha', 'beta'], index=xind)
####
init_profile=init_profile.iloc[::50,:]
init_profile.reset_index(drop=True, inplace=True)
####

column = 'values'
para = pd.Series(name = column, dtype=np.float64)
para = para.astype('object')

# System-level 
para.at['problem'] = 'HeatConduction'
para.at['SpatialDiscretize'] = 'CenteredDifferencing'
para.at['TimeDiscretize'] = 'BackwardEular'
para.at['ODEsolver'] = 'NewtonIteration'
para.at['linearSolver'] = 'numpy linalg'
para.at['CPU'] = 1
para.at['NNmodel']= torch.load('C:/Users/aleks/Projects/ML-student-projects/PyTorch/Model.pt')
para.at['NNmodel'].eval()
para.at['scaling']=pd.read_csv('C:/Users/aleks/Projects/ML-student-projects/PyTorch/data_scaling.csv'\
                                 , index_col=(0))
# Material
para.at['material'] = 'steel'
para.at['material function'] = 'constant'
para.at['density'] = 7850
para.at['conductivity'] = 60.5
para.at['heatCapacity'] = 434

# Grid
para.at['length'] = 1
para.at['numberOfNode'] = len(init_profile)
para.at['x']=init_profile['x']
para.at['scaling']=pd.read_csv('C:/Users/aleks/Projects/ML-student-projects/PyTorch/data_scaling.csv'\
                             , index_col=(0))
# Solution
para.at['numberOfTimeStep'] = 40#400
para.at['deltaTime'] = 0.2
para.at['maxIteration'] = 20
para.at['convergence'] = 1e-2#1E-10
para.at['relaxation'] = 1 # value in [0-1] Very sensitive!!!

# Initial conditions
para.at['InitTeProfile'] = init_profile['Te']
para.at['InitneProfile'] = init_profile['ne']
para.at['InitgradTeProfile'] = init_profile['gradTe']
para.at['InitZbarProfile'] = init_profile['Zbar']
para.at['InitKnProfile'] = init_profile['Kn']


# Boundary conditions
para.at['x=0 type'] = 'heatFlux' #'heatFlux' or 'fixedTemperature'
para.at['x=0 value'] = 0
para.at['x=L type'] = 'heatFlux' #'heatFlux' or 'fixedTemperature'
para.at['x=L value'] = 0  




rho = para['density']
hcp = para['heatCapacity']

dt = para['deltaTime']
length = para['length']
numberOfNode = para['numberOfNode']



dx = length / (numberOfNode -1)

# BC informations
typeX0 = para['x=0 type']
valueX0 = para['x=0 value']
typeXL = para['x=L type']
valueXL = para['x=L value']



####
T=para['InitTeProfile']
x=para['x']
Z=para['InitZbarProfile']
ne=para['InitneProfile']
Kn=para['InitKnProfile']
scale=para['scaling']
k=(Z+0.24)/(Z*(Z+4.2))
gradT=np.gradient(T,x)

#Parameters given by NN:
    #number of points in one domain / number of features (T, gradt, Z, n, kn, x)
Tscaled=(T-scale['T'].loc['mean'])/scale['T'].loc['std']
gradTscaled=(gradT-scale['gradT'].loc['mean'])/scale['gradT'].loc['std']
Zscaled=(Z-scale['Z'].loc['mean'])/scale['Z'].loc['std']
nescaled=(ne-scale['n'].loc['mean'])/scale['n'].loc['std']
Knscaled=(Kn-scale['Kn'].loc['mean'])/scale['Kn'].loc['std']
    

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
alpha, beta = para['NNmodel'](torch.tensor(datavector.values.flatten('F'), dtype=torch.float))
params = pd.DataFrame([[alpha.detach().numpy(),beta.detach().numpy()]],columns=['alpha', 'beta'], index=xind)

for ind,_ in enumerate(x, start=1):
    
    if ind+lng>=len(x)+1:
        break    
    datavector=datavector.drop([0]).reset_index(drop=True)
    datavector.loc[len(datavector)+1]=[Tscaled[ind], gradTscaled[ind],\
                                       Zscaled[ind], nescaled[ind], Knscaled[ind],0]
    xind=np.append(xind, int((lng-1)/2+ind)) #index of the alphabeta
    alpha, beta = para['NNmodel'](torch.tensor(datavector.values.flatten('F'), dtype=torch.float))
    params.loc[xind[-1]]=[alpha.detach().numpy(),beta.detach().numpy()] #??? why indexes are converting to float ???
    
    if ind%5000==0:
            print(f"alpha and beta are calculated for {ind}/{len(x)-lng+1} points") 




#Add alphas and betas at the beginning and end of intervals 
#in order to all arrays  have the same length

for i in range(xind[0]):              
    params = pd.concat([params.iloc[0].to_frame().T.set_index(pd.Index([xind[0]-i-1])), params])
for i in range(len(x)-xind[-1]-1):              
    params = pd.concat([params,params.iloc[0].to_frame().T.set_index(pd.Index([xind[-1]+i+1]))])
    if i%500==0:
            print(f"alpha and beta are calculated for {i}/{len(x)-lng+1} points") 