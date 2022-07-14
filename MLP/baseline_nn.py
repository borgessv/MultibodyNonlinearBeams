# -*- coding: utf-8 -*-
"""
Multilayer Perceptron neural network
"""

import torch
#import numpy as np

class MLP_NN(torch.nn.Module):
  '''Just a salt-of-the-earth MLP'''
  def __init__(self, input_dim, hidden_dim, output_dim, nonlinearity='tanh'):
    super(MLP_NN, self).__init__()
    self.linear1 = torch.nn.Linear(input_dim, hidden_dim)
    self.linear2 = torch.nn.Linear(hidden_dim, hidden_dim)
    self.linear3 = torch.nn.Linear(hidden_dim, output_dim, bias=None)

    for l in [self.linear1, self.linear2, self.linear3]:
      torch.nn.init.orthogonal_(l.weight) # use a principled initialization

    self.nonlinearity = nonlinearity

  def forward(self, x, separate_fields=False):
    h = self.nonlinearity(self.linear1(x))
    h = self.nonlinearity(self.linear2(h))
    #print(self.linear3(h))
    return self.linear3(h)
    


class NN(torch.nn.Module):
    '''Learn arbitrary vector fields that are sums of conservative and solenoidal fields'''
    def __init__(self, input_dim, differentiable_model):
        super(NN, self).__init__()
        self.differentiable_model = differentiable_model
        
    def forward(self, x):
        # traditional forward pass
        return self.differentiable_model(x)

        y = self.differentiable_model(x)
        assert y.dim() == 2 and y.shape[1] == 2, "Output tensor should have shape [batch_size, 2]"
        return y.split(1,1)
    
    def time_derivative(self, x, t=None, separate_fields=False):
        '''NEURAL ODE-STLE VECTOR FIELD'''
        y1 = self.differentiable_model(x)
        #print(y1)
        return y1
        