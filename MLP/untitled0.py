# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 19:36:49 2022

@author: Vitor
"""

from torchdyn.core import NeuralODE
from torchdyn.nn import DataControl, DepthCat, Augmenter, GalLinear, Fourier
from torchdyn.datasets import *
from torchdyn.utils import *

#%load_ext autoreload
#%autoreload 2

d = ToyDataset()
X, yn = d.generate(n_samples=512, dataset_type='moons', noise=.1)

import matplotlib.pyplot as plt

colors = ['orange', 'blue'] 
fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot(111)
for i in range(len(X)):
    ax.scatter(X[i,0], X[i,1], s=1, color=colors[yn[i].int()])