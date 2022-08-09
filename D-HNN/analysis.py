# -*- coding: utf-8 -*-
"""
"""
import os, sys
import torch
import numpy as np
import scipy.integrate
import scipy.io as sio
import matplotli.pyplot as plt
import matlab.engine

from models import MLP
from train import train, get_args
from data_maker import make_dataset


def model_args():
  model_param = {'model': 'FOM',
                 'DoF': {'OutBend'},
                 'gravity': 'GravityOn',
                 'beam_data': 'beam_data_test.xlsx',
                 'tspan': [0,10],
                 'timestep': 0.1,
                 'p0span': [0,0],
                 'q0span': [-np.pi/3,np.pi/3]}
  return model_param

# Dataset:
model_param = model_args()
data, data_raw = make_dataset(model_param, samples=50, train_split=0.5)
n_DoF = data['n_DoF']
M = data['M']
I = data['I']
K = data['K']
C = data['C']

#%% Train NN model:
args = get_args()
args.input_dim = data['x'].shape[1] + 1
args.output_dim = data['dx'].shape[1]
mlp_model = MLP(args.input_dim, args.output_dim, args.hidden_dim)
mlp_results = train(mlp_model, args, data)

# Run simulation using NN model:
def integrate_model(model, t_span, y0, use_torch=True, **kwargs):
    
    def fun(t, x):
        if use_torch:
            x = torch.tensor(x, requires_grad=True, dtype=torch.float32).reshape(1,2*n_DoF)
            t = torch.zeros_like(x[...,:1])
        else:
            x = x.reshape(1,2*n_DoF)
            t = np.zeros_like(x[...,:1])
        dx = model(x, t=t).reshape(-1)
        if use_torch:
            dx = dx.data.numpy()
        return dx

    return scipy.integrate.solve_ivp(fun=fun, t_span=t_span, y0=y0, **kwargs)

# MATLAB interface:
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir) 
eng = matlab.engine.start_matlab()
eng.addpath(parent_dir + "\\SimulationFramework", nargout= 0)
eng.addpath(parent_dir + "\\SimulationFramework\\background", nargout= 0)
eng.addpath(parent_dir + "\\SimulationFramework\\background\\utils", nargout= 0)
eng.addpath(parent_dir + "\\SimulationFramework\\background\\CrossSectionData", nargout= 0)

# Test Parameters:
y0 = np.array([0.0,0.0,0.0,0.0,0.0,0.0]) 
# p0 = np.asarray([0])
# q0 = np.asarray([0])
# p0 = (p0 - np.min(data_raw['x'][:,:n_DoF]))/np.ptp(data_raw['x'][:,:n_DoF])
# # q0 = (q0 - np.min(data_raw['x'][:,n_DoF:]))/np.ptp(data_raw['x'][:,n_DoF:])
# y0 = np.array(np.concatenate([p0,q0])).T
t_span = [0,20]
dt = 0.01
tvec = np.array(np.linspace(t_span[0],t_span[1],int((t_span[1]-t_span[0])/dt)+1))
kwargs = {'t_eval': tvec, 'rtol': 1e-10}
X_nn = integrate_model(mlp_model, t_span, y0, **kwargs)

X_nn = X_nn['y'].T
# p_real = (x_sol[:,:n_DoF])*np.ptp(data_raw['x'][:,:n_DoF]) + np.min(data_raw['x'][:,:n_DoF])
# q_real = (x_sol[:,n_DoF:])*np.ptp(data_raw['x'][:,n_DoF:]) + np.min(data_raw['x'][:,n_DoF:])
# q_real = (x_sol[:,n_DoF:])
# xs_sol = np.concatenate((p_real,q_real),axis=1)
sio.savemat('test.mat', {'X_MLP': X_nn})

X_true = eng.simulation(model_param['model'],model_param['DoF'],model_param['gravity'],matlab.double(tvec),M,I,K,C,matlab.double(y0),nargout=2)[0]
X_true = np.array(X_true)

#%% Plot results:
p_nn, q_nn = X_nn[:,:n_DoF], X_nn[:,n_DoF:]
p_true, q_true = X_true[:,:n_DoF], X_true[:,n_DoF:]
dof = 1
plt.plot(q_true[:,dof],p_true[:,dof],'k', q_nn[:,dof],p_nn[:,dof],'r')
plt.xlabel('q [rad]')
plt.ylabel('p [kg*rad/s]')
plt.grid(True)
plt.legend(['Ground True','NN'])


