# -*- coding: utf-8 -*-
"""
Description:
Dataset creation using simulation outputs of the dynamical system modeled by  
the lumped-parameter method combined with the Hamiltonian mechanics framework.
For more details check MATLAB's 'simulation.m' code.

Author: Vitor Borges Santos - borgessv93@gmail.com
Version: 10-July-2022
"""

import os, sys
import numpy as np
import pandas as pd
import matlab.engine

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

# Paths where the MATLAB codes are located: 
eng = matlab.engine.start_matlab()
eng.addpath(parent_dir + "\\SimulationFramework", nargout= 0)
eng.addpath(parent_dir + "\\SimulationFramework\\background", nargout= 0)
eng.addpath(parent_dir + "\\SimulationFramework\\background\\utils", nargout= 0)
eng.addpath(parent_dir + "\\SimulationFramework\\background\\CrossSectionData", nargout= 0)

def make_dataset(model_param, seed=0, samples=20, train_split=0.5, **kwargs):
    print('\rCreating dataset... {:.1f}% done'.format(0), end='')
    
    # Model parameters:    
    model = model_param['model']
    DoF = model_param['DoF']
    gravity = model_param['gravity']
    
    # Initializing structure properties:
    beam_data = model_param['beam_data']
    beam = pd.read_excel(parent_dir + "\\SimulationFramework\\" + beam_data)
    n_DoF = len(DoF)*int(np.sum(np.array(beam.iloc[3,1:])))
    M, I, K, C = eng.structure_properties(beam_data,DoF,nargout=4)
    #K = eng.times(1e-3,K)
    #C = eng.times(1e-1,C)
    #K = matlab.double(K)
    #K[0][0] = 0
    
    # Simulation timestep:
    dt = model_param['timestep']
    t0 = model_param['tspan'][0]
    t1 = model_param['tspan'][1]
    tvec = np.array([np.linspace(t0,t1,int((t1-t0)/dt)+1)])
    teval = matlab.double(tvec)
    
    # Initial condition parameters:
    np.random.seed(seed)
    p0span = model_param['p0span']
    q0span = model_param['q0span']
    
    # Dataset creation loop:
    X_data, Xdot_data, t_data = [], [], []
    for i in range(samples):
        p0 = (p0span[1]-p0span[0])*np.random.rand(n_DoF) + p0span[0] 
        q0 = (q0span[1]-q0span[0])*np.random.rand(n_DoF) + q0span[0]
        X0 = np.array([np.concatenate([p0,q0])]).T
        X0 = matlab.double(X0)
        
        X, Xdot = eng.simulation(model,DoF,gravity,teval,M,I,K,C,X0,nargout=2)
        
        X = np.array(X)        
        X_data.append(X)
        
        Xdot = np.array(Xdot)        
        Xdot_data.append(Xdot)
             
        t_data.append(tvec)
        
        progress_msg = '\rCreating dataset... {:.1f}% done'.format(100*(i+1)/samples)
        print(progress_msg + '\n' if i == samples-1 else progress_msg, end='')
    x_aux = np.concatenate(X_data)
    t = np.zeros_like(x_aux[...,:1])
    data_raw = {'x': np.concatenate(X_data), 
                'dx': np.concatenate(Xdot_data), 
                't': t}
    
    # Normalizing dataset [-1,1]:
    data = {'x': np.concatenate(X_data), 
            'dx': np.concatenate(Xdot_data), 
            't': t}
    # data['x'][:,:n_DoF] = (data['x'][:,:n_DoF] - np.min(data['x'][:,:n_DoF]))/np.ptp(data['x'][:,:n_DoF])    
    # data['x'][:,n_DoF:] = (data['x'][:,n_DoF:] - np.min(data['x'][:,n_DoF:]))/np.ptp(data['x'][:,n_DoF:]) 
    # data['dx'][:,:n_DoF] = (data['dx'][:,:n_DoF] - np.min(data['dx'][:,:n_DoF]))/np.ptp(data['dx'][:,:n_DoF])    
    # data['dx'][:,n_DoF:] = (data['dx'][:,n_DoF:] - np.min(data['dx'][:,n_DoF:]))/np.ptp(data['dx'][:,n_DoF:])
    
    # make a train/test split
    split_ix = int(len(data['x'])*train_split)
    split_data = {}
    for k in ['x', 'dx','t']:
        split_data[k], split_data[k + '_test'] = data[k][:split_ix], data[k][split_ix:]
    data = split_data
    data.update({'n_DoF': n_DoF, 'M': M, 'I': I, 'K': K, 'C': C})
    
    return data, data_raw
#data,data_raw = make_dataset(samples=2)