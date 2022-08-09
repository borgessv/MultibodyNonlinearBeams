# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 14:07:03 2022

@author: Vitor
"""

from torchdyn.core import NeuralODE
from torchdyn.datasets import *
from torchdyn import *

import torch
import torch.utils.data as data
device = torch.device("cpu")

from data_maker import make_dataset


dataset, dataset_raw = make_dataset(samples=30)

X = dataset_raw['x']
yn = dataset_raw['dx']
X_train = torch.Tensor(X).to(device)
y_train = torch.Tensor(yn).to(device)
train = data.TensorDataset(X_train, y_train)
trainloader = data.DataLoader(train, batch_size=len(X), shuffle=True)

import torch.nn as nn
import pytorch_lightning as pl

class Learner(pl.LightningModule):
    def __init__(self, t_span:torch.Tensor, model:nn.Module):
        super().__init__()
        self.model, self.t_span = model, t_span

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        x, y = batch
        t_eval, y_hat = self.model(x, t_span)
        y_hat = y_hat[-1] # select last point of solution trajectory
        loss = nn.CrossEntropyLoss()(y_hat, y)
        return {'loss': loss}

    def configure_optimizers(self):
        return torch.optim.Adam(self.model.parameters(), lr=0.01)

    def train_dataloader(self):
        return trainloader
    
f = nn.Sequential(
        nn.Linear(20, 400),
        nn.Tanh(),
        nn.Linear(400, 20)
    )
t_span = torch.linspace(0, 10, 100)

model = NeuralODE(f, sensitivity='adjoint', solver='dopri5').to(device)

learn = Learner(t_span, model)
trainer = pl.Trainer(min_epochs=200, max_epochs=300)
trainer.fit(learn)

t_span = torch.linspace(0,10,100)
t_eval, trajectory = model(X_train, t_span)
trajectory = trajectory.detach().cpu()