# -*- coding: utf-8 -*-
"""
"""
import time, copy, os, sys
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matlab.engine

from models import MLP, HNN, DHNN
from train import train
from data_maker import make_dataset
from utils import matlab_interface, ObjectView, normalize, standardize, integrate_model
matlab_eng = matlab_interface()

# plots are shown in the plot pane:
from IPython import get_ipython
ipython = get_ipython()

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)
if 'windows' in sys.platform:
    bar = '\\'
else:
    bar = '/'

#%% Model parameters and dataset creation:
def model_args(as_dict=False):
  model_dict = {'matlab_engine': matlab_eng,
                'model': 'FOM',
                'DoF': {'OutBend'},
                'gravity': 'GravityOn',
                'beam_data': 'beam_data_test.xlsx',
                'disp_progress': 'False',
                'tspan': [0,10],
                'timestep': 0.01,
                'p0span': [0,0],
                'q0span': [-np.pi/9,np.pi/9],
                'train_split': 0.7,
                'normalize': False,
                'norm_range': [0,1],
                'standardize': False}
  return model_dict if as_dict else ObjectView(model_dict)

# Dataset:
model_param = model_args()
data_raw = make_dataset(model_param, samples=50, train_split=0.5)
n_DoF, M, I, K, C = data_raw['n_DoF'], data_raw['M'], data_raw['I'], data_raw['K'], data_raw['C']

#%% Preprocessing:
if model_param.normalize:
    data = normalize(copy.deepcopy(data_raw), model_param.norm_range)
elif model_param.standardize:
    data = standardize(copy.deepcopy(data_raw))
else:
    data = copy.deepcopy(data_raw)

# Train/test split:
split_ix = int(len(data['x'])*model_param.train_split)
split_data = {}
for k in ['x', 'dx','t']:
    split_data[k], split_data[k + '_test'] = data[k][:split_ix], data[k][split_ix:]
data = split_data

#%% Train NN model:
if 'dhnn_model' and 'mlp_model' and 'hnn_model' in locals():
    del dhnn_model, mlp_model, hnn_model
    
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

# hnn_model = HNN(args.input_dim, args.hidden_dim)
# hnn_results = train(hnn_model, args, data)

dhnn_model = DHNN(args.input_dim, args.hidden_dim)
dhnn_results = train(dhnn_model, args, data)

#%% Test simulation model:
print('\rRunning test simulations... {:.1f}% done'.format(0), end='')

# Test Parameters:
tspan = [0,15]
dt = 0.01
tvec = np.array(np.linspace(tspan[0],tspan[1],int((tspan[1]-tspan[0])/dt)+1))
kwargs = {'t_eval': tvec, 'rtol': 1e-10}
p0span = model_param.p0span
q0span = model_param.q0span
n_test = 1

tsim_mlp, tsim_dhnn, tsim_true, t_mlp_true, t_dhnn_true = [], [], [], [], []
for i in range(n_test):
    p0_test = 0*((p0span[1]-p0span[0])*np.random.rand(n_DoF) + p0span[0] )
    q0_test = 0*((q0span[1]-q0span[0])*np.random.rand(n_DoF) + q0span[0])
    X0_test_raw = np.concatenate([p0_test,q0_test])

    if model_param.normalize:
        norm_range = model_param.norm_range
        p0_test = np.ptp(norm_range)*(p0_test - np.min(data_raw['x'][:,:n_DoF]))/np.ptp(data_raw['x'][:,:n_DoF]) + norm_range[0]   
        q0_test = np.ptp(norm_range)*(q0_test - np.min(data_raw['x'][:,n_DoF:]))/np.ptp(data_raw['x'][:,n_DoF:]) + norm_range[0]
        X0_test = np.concatenate([p0_test,q0_test])
    elif model_param.standardize:
        p0_test = (p0_test - np.mean(data_raw['x'][:,:n_DoF]))/np.std(data_raw['x'][:,:n_DoF])   
        q0_test = (q0_test - np.mean(data_raw['x'][:,n_DoF:]))/np.std(data_raw['x'][:,n_DoF:])
        X0_test = np.concatenate([p0_test,q0_test])
    else:
        X0_test = X0_test_raw   

    # Run MLP model simulation for the test case
    t = time.time()
    X_mlp = integrate_model(mlp_model, tspan, X0_test, n_DoF, **kwargs) # MLP NN model simulation
    tsim_mlp.append(time.time() - t)

    # Run DHNN model simulation for the test case
    t = time.time()
    X_dhnn = integrate_model(dhnn_model, tspan, X0_test, n_DoF, **kwargs) # DHNN model simulation
    tsim_dhnn.append(time.time() - t)
 
    # Run MATLAB simulation for the test case:
    t = time.time()
    X_true = matlab_eng.simulation(model_param.model,model_param.DoF,model_param.gravity,
                                   matlab.double(tvec),matlab.double(X0_test_raw),
                                   model_param.disp_progress,nargout=2)[0]
    tsim_true.append(time.time() - t)
    
    t_mlp_true.append((tsim_mlp[i]/tsim_true[i])*100)
    t_dhnn_true.append((tsim_dhnn[i]/tsim_true[i])*100)
    
    progress_msg = '\rrunning test simulations... {:.1f}% done'.format(100*(i+1)/n_test)
    print(progress_msg + '\n' if i == n_test-1 else progress_msg, end='')
    
# Denormalize results (if needed):
if model_param.normalize:
    X_mlp_norm = X_mlp['y'].T
    p_mlp = (X_mlp_norm[:,:n_DoF] - norm_range[0])*np.ptp(data_raw['x'][:,:n_DoF])/np.ptp(norm_range) + np.min(data_raw['x'][:,:n_DoF])
    q_mlp = (X_mlp_norm[:,n_DoF:] - norm_range[0])*np.ptp(data_raw['x'][:,n_DoF:])/np.ptp(norm_range) + np.min(data_raw['x'][:,n_DoF:])
    X_MLP = np.array(np.concatenate([p_mlp.T,q_mlp.T])).T
    X_dhnn_norm = X_dhnn['y'].T
    p_dhnn = (X_dhnn_norm[:,:n_DoF] - norm_range[0])*np.ptp(data_raw['x'][:,:n_DoF])/np.ptp(norm_range) + np.min(data_raw['x'][:,:n_DoF])
    q_dhnn = (X_dhnn_norm[:,n_DoF:] - norm_range[0])*np.ptp(data_raw['x'][:,n_DoF:])/np.ptp(norm_range) + np.min(data_raw['x'][:,n_DoF:])
    X_DHNN = np.array(np.concatenate([p_dhnn.T,q_dhnn.T])).T
elif model_param.standardize:
    X_mlp_std = X_mlp['y'].T
    p_mlp = X_mlp_std[:,:n_DoF]*np.std(data_raw['x'][:,:n_DoF]) + np.mean(data_raw['x'][:,:n_DoF])
    q_mlp = X_mlp_std[:,n_DoF:]*np.std(data_raw['x'][:,n_DoF:]) + np.mean(data_raw['x'][:,n_DoF:])
    X_MLP = np.array(np.concatenate([p_mlp.T,q_mlp.T])).T
    X_dhnn_std = X_dhnn['y'].T
    p_dhnn = X_dhnn_std[:,:n_DoF]*np.std(data_raw['x'][:,:n_DoF]) + np.mean(data_raw['x'][:,:n_DoF])
    q_dhnn = X_dhnn_std[:,n_DoF:]*np.std(data_raw['x'][:,n_DoF:]) + np.mean(data_raw['x'][:,n_DoF:])
    X_DHNN = np.array(np.concatenate([p_dhnn.T,q_dhnn.T])).T
else:
    X_MLP = X_mlp['y'].T
    X_DHNN = X_dhnn['y'].T

X_true = np.array(X_true)

# Save NN simulation results:
sio.savemat(parent_dir + bar + "SimulationFramework" + bar + "test_results" + bar + "test.mat", {'X_FOM': X_true, 'X_MLP': X_MLP, 'X_DHNN': X_DHNN, 'tvec': matlab.double(tvec)})

# NN relative computaional cost:
print('\nComputational cost relative to MATLAB simulation:\nMLP: {:.1f}% +/- {:.2f}% \nDHNN: {:.1f}% +/- {:.2f}%'.format(np.mean(t_mlp_true),np.std(t_mlp_true),np.mean(t_dhnn_true),np.std(t_dhnn_true)))


#%% Phase portrait for last test run:
ipython.magic('matplotlib inline')
dof = 0 # DoF considered for the phase portrait 
p_mlp, q_mlp = np.split(X_MLP, 2, axis=1)
p_dhnn, q_dhnn = np.split(X_DHNN, 2, axis=1)
p_true, q_true = np.split(X_true, 2, axis=1)

fig = plt.figure()
plt.plot(q_true[:,dof],p_true[:,dof],'k', q_mlp[:,dof],p_mlp[:,dof],'--b', q_dhnn[:,dof],p_dhnn[:,dof],'--r')
plt.xlabel('q [rad]')
plt.ylabel('p [kg*rad/s]')
plt.grid(True)
plt.legend(['Ground Truth','MLP','DHNN'],loc='upper right')
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