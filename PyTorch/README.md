# This folder contains the key files for solving the non-local transport problem
## Generate the data neseccary for neural network to train
We used local data that has been generate via `python3 generate_QimpactTrainingData.py 100 0.01 0.2 40000 105`
However, I forgot to include these data to `gitignore.git` and have accidentally uploaded them to this repository
`NN-heatflux.ipynb` and `NN_training.py` contain training part of the neural network, whereas `HeatFluxModel.py` and `HeatFluxData.py` contain the definiotion of the NN itself and data preparation respectevely

# Time evolution of the nonlocal problem 
The script responsible for time evolution problem can be reached either via 'time-evolv-heatflux.ipynb' notebook, or can be started directly by running `parameter.py` file

All parameters of the problem, such as number of used profiles points, boundary conditions, definition of parameters $\alpha$ and $\beta$, number of timesteps and its magnitude $dt$ are contained and may be modified in `parameter.py` file.

`heatConduction.py` package contains necessary components to solve the PDEs such as calculations of jacobians, Newton iterations, second order derivatives, etc.
`postprocessing.py` package is responsible for visualization of the results.

# PyTorch tutorials for NTH-ML


1. PyTorch Lightning Tutorial, https://becominghuman.ai/pytorch-lightning-tutorial-1-getting-started-5f82e06503f6
2. Introduction to PyTorch Lightning, https://pytorch-lightning.readthedocs.io/en/stable/
3. Regression using PyTorch Lightning, "Bike Share Regression PyTorch Lightning.ipynb", https://github.com/shotleft/how-to-python.git

# !!!
The structure of the modules solving the time evolution of temperature distribution (`heatConduction.py`, `postprocessing.py`,  `parameter.py`)  are based on github heatConduction repository published by https://github.com/rickfu415
