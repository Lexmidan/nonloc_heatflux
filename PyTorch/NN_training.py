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


def train_model(numEpochs):
    training_data = './Data/scaled_QdataKn23width100microns.csv'
    data_scaling = pd.read_csv('./Data/data_scaling.csv', header=0, index_col=0)
    test_split = 0.5
    train_split = 0.8
    dropouts = ['Kn', 'x']
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
    trainer = pl.Trainer(max_epochs = numEpochs)
    trainer.fit(model, train_loader, validation_loader)
    
    return model