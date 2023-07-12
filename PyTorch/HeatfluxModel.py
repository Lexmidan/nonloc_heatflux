import torch
import pytorch_lightning as pl
from pytorch_lightning import Trainer
import numpy as np
import parameter as para
import heatConduction as hc
import scipy.constants as const
import scipy

class HeatFluxModel(pl.LightningModule): 

### Model ###
    def __init__(self, Nfeatures, Nlayers, Ntargets, scaling, Nfields):
        super(HeatFluxModel, self).__init__() # TODO: if not "cannot assign module before Module.__init__() call"
        # Initialize layers
        self.fcIn = torch.nn.Linear(Nfeatures, Nlayers[0])
        self.fc = []
        for i in range(len(Nlayers)-1):
            self.fc.append(torch.nn.Linear(Nlayers[i], Nlayers[i+1]))
        self.fcOut = torch.nn.Linear(Nlayers[len(Nlayers)-1], Ntargets)
        # Keep data scaling
        self.scaling = scaling
        # TODO: improve evaluation of Npoints
        self.Npoints = int(Nfeatures / Nfields)
        # TODO: better place to define mse_loss
        self.mse_loss = torch.nn.MSELoss(reduction = 'mean')

    def forward(self, x):
        x = torch.relu(self.fcIn(x))
        for i in range(len(self.fc)):
            x = torch.relu(self.fc[i](x))
        x = self.fcOut(x)
        return x

### The Optimizer ### 
    def configure_optimizers(self):
        #optimizer = torch.optim.Adam(self.parameters(), lr=0.05)#l_rate) # TODO: should be a parameter
        optimizer = torch.optim.SGD(self.parameters(), lr=0.05)#l_rate) # TODO: should be a parameter
        return optimizer

### Training ### 
    def training_step(self, batch, batch_idx):
        x, y = batch
        # Evaluate physical model using data scaling
        logits = self.forward(x)
        # Evaluate loss comparing to the kinetic heat flux in y
        loss = self.mse_loss(logits, y)
        # Add logging
        logs = {'loss': loss}
        return {'loss': loss, 'log': logs}

### Validation ### 
    def validation_step(self, batch, batch_idx):
        x, y = batch
        # Evaluate physical model using data scaling
        logits = self.forward(x)
        # Evaluate loss comparing to the kinetic heat flux in y
        loss = self.mse_loss(logits, y)
        return {'val_loss': loss}

    # Define validation epoch end
    def validation_epoch_end(self, outputs):
        avg_loss = torch.stack([x['val_loss'] for x in outputs]).mean()
        tensorboard_logs = {'val_loss': avg_loss}
        return {'avg_val_loss': avg_loss, 'log': tensorboard_logs}


class AlphaBetaModel(HeatFluxModel): 

    def __init__(self, Nfeatures, Nlayers, scaling, Nfields): # AlphaModel uses alpha, beta outputs (2)
        super(AlphaBetaModel, self).__init__(Nfeatures, Nlayers, 2, scaling, Nfields)

    def scaled_heatflux_model(self, x):
        Q = self.forward(x)[:, 0]
        return Q

    def scaled_beta_model(self, x):
        nln = self.forward(x)[:, 1]
        return nln
    
    def heatflux_model(self, x):
        mean = self.scaling['Q']['mean']
        std = self.scaling['Q']['std']
        return self.scaled_heatflux_model(x) * std + mean    

    def beta_model(self, x):
        # TODO: beta data should also have scaling, but it does not yet
        mean = self.scaling['beta']['mean']
        std = self.scaling['beta']['std']
        beta = self.scaled_beta_model(x) * std + mean
        # Make sure the power of diffusivity is positive
        beta[beta<1e-6] = 1e-6
        return beta    
    
    def alpha_model(self, x):
        feature = self.feature_xc(x)
        Z = feature['Z']; T =  feature['T']; gradT = feature['gradT'], ne = feature['n'], 
        m_e = 9.1094*1e-28
        q_e = 4.8032*1e-10 
        coulog = 23-np.log(np.sqrt(ne)*Z/T**1.5)
        Kb=1.602178e-12
        v=np.sqrt(T*Kb/m_e)
        Gamma= 4 * const.pi * q_e**4/m_e**2
        lamb = v**4/(ne*Gamma*coulog)*1/np.sqrt(Z+1)
        Kn = -lamb*gradT/T
        #Kn_nonloc=scipy.signal.convolve(Kn, hc.gaussian_kernel(size = self.fcIn.in_features/4, sigma = 6), mode='same')
        #alpha = hc.calc_alpha(self.heatflux_model(x), self.beta_model(x), Z, T, gradT, Kn_nonloc) 
        diffusive_heatflux_model = self.local_heatflux_beta_model(x)
        diffusive_heatflux_model[diffusive_heatflux_model<1e-3]=self.heatflux_model(x)[diffusive_heatflux_model<1e-3]
        alpha = self.heatflux_model(x) / diffusive_heatflux_model
        alpha[alpha<1e-6] = 1e-6
        alpha = hc.alpha_cor(alpha,Kn)
        
        # #diffusive_heatflux_model = self.local_heatflux_model(x)
        # diffusive_heatflux_model = self.local_heatflux_beta_model(x)
        # #diffusive_heatflux_model[diffusive_heatflux_model<1e-4]=1e11
        # alpha = self.heatflux_model(x) / diffusive_heatflux_model
        # #alpha[diffusive_heatflux_model<1e-4]=10
        return alpha

    def local_heatflux_beta_model(self, x):
        # Extract features for modified local heat flux
        feature = self.feature_xc(x)
        # Nonlinear conductivity power
        beta = self.beta_model(x)
        # Get beta local heat flux
        Z = feature['Z']; TkeV = 1e-3 * feature['T']; gradTkeV = 1e-3 * feature['gradT']
        kQSH = 6.1e+02 * 1e3**2.5 * 1e3 # scaling constant consistent with SCHICK and T in keV
        ### TODO
        # Workaround to get proper tensor dimensions
        local_heatflux_beta_model = - kQSH / Z[:] * ((Z[:] + 0.24) / (Z[:] + 4.2))\
          * TkeV[:]**beta[:] * gradTkeV[:]
        ###        
        return local_heatflux_beta_model
    
    def local_heatflux_model(self, x):
        # Extract features for modified local heat flux
        # The order MUST correspond to 'generate_QimpactTrainingData.py' T, gradT, Z, n
        featureName = ['T', 'gradT', 'Z']
        feature = {}
        ip = int(self.Npoints / 2)
        i = 0
        for name in featureName:
            mean, std = self.scaling[name]['mean'], self.scaling[name]['std']
            # Batch values of the feature at xc
            feature[name] = x[:, i * self.Npoints + ip] * std + mean
            i = i+1        
        # Get local heat flux
        Z = feature['Z']; T = feature['T']; gradT = feature['gradT']
        kQSH = 6.1e+02 # scaling constant consistent with SCHICK and T in keV
        ### TODO
        # Workaround to get proper tensor dimensions
        local_heatflux_model = - kQSH / Z[:] * ((Z[:] + 0.24) / (Z[:] + 4.2))\
          * T[:]**2.5 * gradT[:]
        ###
        return local_heatflux_model
    
    # TODO: make this more independent from the dataset order
    def feature_xc(self, x):
         # The order MUST correspond to 'generate_QimpactTrainingData.py' T, gradT, Z, n
        featureName = ['T', 'gradT', 'Z', 'n']
        feature = {}
        ip = int(self.Npoints / 2)
        i = 0
        for name in featureName:
            mean, std = self.scaling[name]['mean'], self.scaling[name]['std']
            # Batch values of the feature at xc
            feature[name] = x[:, i * self.Npoints + ip] * std + mean
            i = i+1
        return feature
