# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 01:24:03 2022

@author: Vitor
"""
import torch, time, sys
import autograd
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
import scipy.integrate
solve_ivp = scipy.integrate.solve_ivp

EXPERIMENT_DIR = './experiment-pend'
sys.path.append(EXPERIMENT_DIR)

#from data import get_dataset, get_field, get_trajectory, dynamics_fn, hamiltonian_fn
from nn_models import MLP
from hnn import HNN
#from utils import L2_loss

RK4 = ''
def get_args():
    return {'input_dim': 2,
         'hidden_dim': 200,
         'learn_rate': 1e-3,
         'nonlinearity': 'relu',
         'total_steps': 2000,
         'field_type': 'solenoidal',
         'print_every': 200,
         'name': 'pend',
         'gridsize': 10,
         'input_noise': 0.5,
         'seed': 0,
         'save_dir': './{}'.format(EXPERIMENT_DIR),
         'fig_dir': './figures'}


class ObjectView(object):   
    def __init__(self, d): self.__dict__ = d
    pass
    
args = ObjectView(get_args())

def get_model(args, baseline):
    output_dim = args.input_dim if baseline else 2
    nn_model = MLP(args.input_dim, args.hidden_dim, output_dim, args.nonlinearity)
    model = HNN(args.input_dim, differentiable_model=nn_model,field_type=args.field_type, baseline=baseline)
    
    model_name = 'baseline' if baseline else 'hnn'
    path = "pend{}-{}.tar".format(RK4, model_name)
    model.load_state_dict(torch.load(path))
    return model

def integrate_model(model, t_span, y0, **kwargs):
    
    def fun(t, np_x):
        x = torch.tensor( np_x, requires_grad=True, dtype=torch.float32).view(1,2)
        dx = model.time_derivative(x).data.numpy().reshape(-1)
        return dx
    
    return solve_ivp(fun=fun, t_span=t_span, y0=y0, **kwargs)

# integrate along those fields starting from point (1,0)
base_model = get_model(args, baseline=True)
dt = 0.1
tf = 10
#tspan = [x*dt for x in range(int(tf/dt + 1))]
t_span = [0,10]
q0max = np.pi
p0max = 0
q0 = q0max*(2*np.random.rand(1)-1) 
p0 = p0max*(2*np.random.rand(1)-1)
y0 = np.concatenate([p0, q0])
kwargs = {'t_eval': np.linspace(t_span[0], t_span[1], 100), 'rtol': 1e-10}
base_ivp = integrate_model(base_model, t_span, y0, **kwargs)