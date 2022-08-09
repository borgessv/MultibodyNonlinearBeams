# -*- coding: utf-8 -*-
"""
"""
import time
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matlab.engine

from models import MLP
from train import train
from data_maker import make_dataset
from utils import *
matlab_eng = matlab_interface()

#%% Model parameters and dataset creation:
def model_args(as_dict=False):
  model_dict = {'matlab_engine': matlab_eng,
                'model': 'FOM',
                'DoF': {'OutBend'},
                'gravity': 'GravityOn',
                'beam_data': 'beam_data_test.xlsx',
                'disp_progress': 'False',
                'tspan': [0,20],
                'timestep': 0.01,
                'p0span': [0,0],
                'q0span': [-np.pi/6,np.pi/6],
                'normalize': False,
                'norm_range': [-1,1],
                'standardize': True}
  return model_dict if as_dict else ObjectView(model_dict)

# Dataset:
model_param = model_args()
data_raw = make_dataset(model_param, samples=50, train_split=0.5)
n_DoF, M, I, K, C = data_raw['n_DoF'], data_raw['M'], data_raw['I'], data_raw['K'], data_raw['C']

# Preprocessing:
if model_param.normalize:
    data = normalize(data_raw, model_param.norm_range)
elif model_param.standardize:
    data = standardize(data_raw)
else:
    data = data_raw

# Train/test split:
split_ix = int(len(data['x'])*train_split)
split_data = {}
for k in ['x', 'dx','t']:
    split_data[k], split_data[k + '_test'] = data[k][:split_ix], data[k][split_ix:]
data = split_data

#%% Train NN model:
def train_args(as_dict=False):
  train_dict = {'input_dim': 2*n_DoF+1,
                'hidden_dim': 50, # capacity
                'output_dim': 2*n_DoF,
                'learning_rate': 1e-3, 
                'test_every': 1,
                'print_every': 200,
                'batch_size': 1001,
                'total_steps': 10000,  # because we have a synthetic dataset
                'device': 'cpu', # {"cpu", "cuda"} for using GPUs
                'seed': 0,
                'as_separate': False,
                'decay': 0}
  return train_dict if as_dict else ObjectView(train_dict)

args = train_args()
mlp_model = MLP(args.input_dim, args.output_dim, args.hidden_dim)
mlp_results = train(mlp_model, args, data)

#%% Test simulation model:
print('\rRunning test simulations... {:.1f}% done'.format(0), end='')

# Test Parameters:
tspan = [0,10]
dt = 0.01
tvec = np.array(np.linspace(tspan[0],tspan[1],int((tspan[1]-tspan[0])/dt)+1))
kwargs = {'t_eval': tvec, 'rtol': 1e-10, 'method': 'LSODA'}
p0span = model_param.p0span
q0span = model_param.q0span
n_test = 1

tsim_nn, tsim_true, tsim_ratio = [], [], []
for i in range(n_test):
    p0_test = data_raw['x'][0,:n_DoF]#(p0span[1]-p0span[0])*np.random.rand(n_DoF) + p0span[0] 
    q0_test = data_raw['x'][0,n_DoF:]#(q0span[1]-q0span[0])*np.random.rand(n_DoF) + q0span[0]
    X0_test_raw = np.concatenate([p0_test,q0_test])
    
    if model_param.normalize:
        norm_range = model_param.norm_range
        p0_test = np.ptp(norm_range)*(p0_test - np.min(data_raw['x'][:,:n_DoF]))/np.ptp(data_raw['x'][:,:n_DoF]) + norm_range[0]   
        q0_test = np.ptp(norm_range)*(q0_test - np.min(data_raw['x'][:,n_DoF:]))/np.ptp(data_raw['x'][:,n_DoF:]) + norm_range[0]
        X0_test = np.concatenate([p0_test,q0_test])
    else:
        X0_test = X0_test_raw   

    # Run NN model simulation for the test case
    t = time.time()
    X = integrate_model(mlp_model, tspan, X0_test, n_DoF, **kwargs) # NN model simulation
    k = X['y'].T
    tsim_nn.append(time.time() - t)

    # Run MATLAB simulation for the test case:
    t = time.time()
    X_true = matlab_eng.simulation(model_param.model,model_param.DoF,model_param.gravity,
                                   matlab.double(tvec),M,I,K,C,matlab.double(X0_test_raw),
                                   model_param.disp_progress,nargout=2)[0]
    tsim_true.append(time.time() - t)
    
    tsim_ratio.append((tsim_nn[i]/tsim_true[i])*100)
    
    progress_msg = '\rRunning test simulations... {:.1f}% done'.format(100*(i+1)/n_test)
    print(progress_msg + '\n' if i == n_test-1 else progress_msg, end='')
    
# Rearrange results from last test run:
    
if model_param.normalize:
    X_nn_norm = X['y'].T
    p_nn = (X_nn_norm[:,:n_DoF] - norm_range[0])*np.ptp(data_raw['x'][:,:n_DoF])/np.ptp(norm_range) + np.min(data_raw['x'][:,:n_DoF])
    q_nn = (X_nn_norm[:,n_DoF:] - norm_range[0])*np.ptp(data_raw['x'][:,n_DoF:])/np.ptp(norm_range) + np.min(data_raw['x'][:,n_DoF:])
    X_nn = np.array(np.concatenate([p_nn.T,q_nn.T])).T
else:
    X_nn = X['y'].T

X_true = np.array(X_true)

# Save NN simulation results:
sio.savemat('test.mat', {'X_MLP': X_nn})

# NN relative computaional cost:
print('NN Relative Simulation Computational Cost: {:.1f}% +/- {:.2f}%'.format(np.mean(tsim_ratio),np.std(tsim_ratio)))

#%% Phase portrait from last test run:
dof = -1 # DoF considered for the phase portrait 
p_nn, q_nn = np.split(X_nn, 2, axis=1)
p_true, q_true = np.split(X_true, 2, axis=1)

fig = plt.figure()
plt.plot(q_true[:,dof],p_true[:,dof],'k', q_nn[:,dof],p_nn[:,dof],'r')
plt.xlabel('q [rad]')
plt.ylabel('p [kg*rad/s]')
plt.grid(True)
plt.legend(['Ground Truth','NN'],loc='upper right')
plt.show()




#y0 = np.array([0.0,0.0,0.0,0.0,0.0,0.0]) 
# p0 = np.asarray([0])
# q0 = np.asarray([0])
# p0 = (p0 - np.min(data_raw['x'][:,:n_DoF]))/np.ptp(data_raw['x'][:,:n_DoF])
# # q0 = (q0 - np.min(data_raw['x'][:,n_DoF:]))/np.ptp(data_raw['x'][:,n_DoF:])
# y0 = np.array(np.concatenate([p0,q0])).T

# p_real = (x_sol[:,:n_DoF])*np.ptp(data_raw['x'][:,:n_DoF]) + np.min(data_raw['x'][:,:n_DoF])
# q_real = (x_sol[:,n_DoF:])*np.ptp(data_raw['x'][:,n_DoF:]) + np.min(data_raw['x'][:,n_DoF:])
# q_real = (x_sol[:,n_DoF:])
# xs_sol = np.concatenate((p_real,q_real),axis=1)