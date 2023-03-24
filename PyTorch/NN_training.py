# -*- coding: utf-8 -*-
"""
NN training part

@author: aleks
"""
import pandas as pd
import HeatfluxModel as hfm
from torch.utils.data import DataLoader
import pytorch_lightning as pl
import HeatfluxData as hfd
import numpy as np
import torch
import HeatfluxData as hfd

def train_model():
    training_data = 'scaled_QdataKn21width100microns.csv'
    data_scaling = pd.read_csv('data_scaling.csv', header=0, index_col=0)
    test_split = 0.5
    train_split = 0.8
    dropouts = []
    scaled_Qdata = pd.read_csv(training_data, header=0, index_col=0)
    Nfields = 6 - len(dropouts) # T, gradT, Z, n, Kn, x
    
    
    _, train_set, validation_set, _ = hfd.heat_flux_datasets(scaled_Qdata, test_split, train_split, dropouts)
    train_loader = DataLoader(dataset = train_set, batch_size = 128)
    validation_loader = DataLoader(dataset = validation_set, batch_size = 128)
    
    # Special object for visualization
    Nfeatures = train_set[0][0].size()[0]# TODO: find a better way than extracting the size via Tensor

    Nlayers = [30, 30]
    NNmodelargs=pd.Series([Nfeatures, Nlayers, data_scaling, Nfields], dtype=object)
    NNmodelargs.to_pickle('./NN/NN_model_args.pkl')
    model = hfm.AlphaBetaModel(*NNmodelargs)
    trainer = pl.Trainer(max_epochs = 100)
    trainer.fit(model, train_loader, validation_loader)
    
    return model

''' TEST
model=trained_model()
init_profile=pd.read_csv('init_profile.csv', index_col=(0))
data_scaling = pd.read_csv(f'data_scaling.csv', header=0, index_col=0)
nescaled  = (init_profile['ne']-data_scaling['n'].loc['mean'])/data_scaling['n'].loc['std']
Zscaled = (init_profile['Zbar']-data_scaling['Z'].loc['mean'])/data_scaling['Z'].loc['std']
Knscaled  = (init_profile['Kn']-data_scaling['Kn'].loc['mean'])/data_scaling['Kn'].loc['std']
Tscaled  = (init_profile['Te']-data_scaling['T'].loc['mean'])/data_scaling['T'].loc['std']
gradTscaled  = (init_profile['gradTe']-data_scaling['gradT'].loc['mean'])/data_scaling['gradT'].loc['std'] 
xdsd=np.array(init_profile['x'])
def get_data_qless(x, T, gradT, Z, n, Kn, width, step):  
    rad=0     #finds out how many points fits in (!)radius(!) range
    s=0
    while s<=width/2:
        s=x[rad]-x[0]
        rad+=1
    
    numPoints = len(Z[0:2*rad:step]) #int(2*rad/step)+1 # the plus one because of integer evaluation
    numFields = 5 #T, gradT, Z, n, Kn
    Qdata=np.empty((0,numFields*numPoints), int) #2 * rad "#of points in interval" * 5 "for each phsy quantity" + 2 "for Q and beta"
    for ind, _ in enumerate(x):  #x_min=x[ind], x_max=x[ind+2*rad], x_c=x[ind+rad]
        datapoint=np.array([])          
        if ind+2*rad>=len(x)+1:
            break    
        else:
            datapoint=np.append(datapoint, T[ind:ind+2*rad:step]) #append all Te in xmin-xmax
            datapoint=np.append(datapoint, gradT[ind:ind+2*rad:step]) #append all gradTe in xmin-xmax
            datapoint=np.append(datapoint, Z[ind:ind+2*rad:step]) #append all Zbar in xmin-xmax
            datapoint=np.append(datapoint, n[ind:ind+2*rad:step]) #append all gradTe in xmin-xmax
            datapoint=np.append(datapoint, Kn[ind:ind+2*rad:step]) #append all Knudsen number in xmin-xmax
            # TODO: what is the appropriate scaling here? Global (max(x)-min(x)) might be to large!
            Qdata=np.append(Qdata,[datapoint], axis=0)
            
            
            if ind%5000==0:
                print(f"We're done with {ind}/{len(x)-2*rad+1} points") 
    #Naming of columns 
    return Qdata, numPoints, rad, (ind-1)+rad #first and last centers of the interval
dataNN, numPoints, fc, lc =get_data_qless(xdsd, Tscaled, gradTscaled, 
                                 Zscaled, nescaled,
                                 Knscaled, 
                                 init_profile['width'][0], int(init_profile['step'][0]))
dataNN =torch.tensor(np.c_[dataNN, np.empty([len(dataNN), numPoints])]).float()
albas=model.alpha_model(dataNN.float()).detach().numpy()
'''