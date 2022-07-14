# -*- coding: utf-8 -*-
"""
Description:
Neural network training framework

Author: Vitor Borges Santos
Version: 11-July-2022

Adapted from: Hamiltonian Neural Networks (2019) by Sam Greydanus, Misko Dzamba
and Jason Yosinski
"""

import torch
import numpy as np
import os, sys
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(PARENT_DIR)

from baseline_nn import MLP_NN, NN
from data_maker import make_dataset
from data import get_dataset

def L2_loss(u, v):
    return (u-v).pow(2).mean()

def train(seed=0, train_steps=2000):
    
    # set random seed
    torch.manual_seed(seed)
    np.random.seed(seed)

    # init model and optimizer 
    n_element = 1
    input_dim = 2*n_element
    hidden_dim = 200
    output_dim = input_dim
    nonlinearity = torch.tanh
    learn_rate = 1e-3
    
    nn_model = MLP_NN(input_dim, hidden_dim, output_dim, nonlinearity)
    model = NN(input_dim, differentiable_model=nn_model)
    optim = torch.optim.Adam(model.parameters(), learn_rate, weight_decay=1e-4)

    # arrange data
    data = make_dataset(seed=seed, samples=20)
    X = torch.tensor(data['X'], requires_grad=True, dtype=torch.float32)
    test_X = torch.tensor(data['test_X'], requires_grad=True, dtype=torch.float32)
    dXdt = torch.Tensor(data['dX'])
    test_dXdt = torch.Tensor(data['test_dX'])

    # vanilla train loop
    print("\rtraining MLP neural network... {:.1f}% done".format(0), end='')
    stats = {'train_loss': [], 'test_loss': []}
    for step in range(train_steps):
        # train step 
        dXdt_hat = model.time_derivative(X)
        loss = L2_loss(dXdt, dXdt_hat)
        loss.backward(); optim.step(); #optim.zero_grad()
    
        # run test data
        test_dXdt_hat = model.time_derivative(test_X)
        test_loss = L2_loss(test_dXdt, test_dXdt_hat)
    
        # logging
        stats['train_loss'].append(loss.item())
        stats['test_loss'].append(test_loss.item())
        
        print("step {}, train_loss {:.4e}, test_loss {:.4e}".format(step, loss.item(), test_loss.item()))

        progress_msg = '\rtraining MLP neural network... {:.1f}% done'.format(100*(step+1)/train_steps)
        #print(progress_msg + '\n' if step == train_steps-1 else progress_msg, end='')
        
    train_dXdt_hat = model.time_derivative(X)
    train_dist = (dXdt - train_dXdt_hat)**2
    test_dXdt_hat = model.time_derivative(test_X)
    test_dist = (test_dXdt - test_dXdt_hat)**2
    print('Final train loss {:.4e} +/- {:.4e}\nFinal test loss {:.4e} +/- {:.4e}'
    .format(train_dist.mean().item(), train_dist.std().item()/np.sqrt(train_dist.shape[0]),
            test_dist.mean().item(), test_dist.std().item()/np.sqrt(test_dist.shape[0])))
  
    return model, stats

if __name__ == "__main__":
    model, stats = train()

    # save
    os.makedirs(THIS_DIR) if not os.path.exists(THIS_DIR) else None
    label = '-MLP' #if args.baseline else '-hnn'
    #label = '-rk4' + label if args.use_rk4 else label
    path = '{}/{}{}.tar'.format(THIS_DIR, 'model', label)
    torch.save(model.state_dict(), path)
