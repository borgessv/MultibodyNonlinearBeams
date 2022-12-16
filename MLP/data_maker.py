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
import scipy.io
import matlab.engine

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

if 'windows' in sys.platform:
    bar = '\\'
else:
    bar = '/'

def make_dataset(model_param, seed=0, samples=20, train_split=0.5, **kwargs):
    print('\rCreating dataset... {:.1f}% done'.format(0), end='')
    
    # Model parameters:  
    eng = model_param.matlab_engine
    model = model_param.model
    DoF = model_param.DoF
    gravity = model_param.gravity
    disp_progress = model_param.disp_progress
    
    # Initializing structure properties:
    beam_data = model_param.beam_data
    beam = pd.read_excel(parent_dir + bar + "SimulationFramework" + bar + beam_data)
    n_DoF = len(DoF)*int(np.sum(np.array(beam.iloc[3,1:])))
    M, I, K, C = eng.structure_properties(beam_data,DoF,disp_progress,nargout=4)
    #x0 = np.array([np.zeros(n_DoF)])
    #Xeq = eng.solve_equilibrium(K,DoF,gravity,'FOM',matlab.double(x0),nargout=1)
    #n_modes = 2
    #phi_r = eng.create_rom(n_modes,DoF,Xeq,nargout=1)
    # K = eng.times(1e-3,K)
    #C = eng.times(1e2,C)
    #K = matlab.double(K)
    #K[0][0] = 0
    #with open(parent_dir + bar + "SimulationFramework" + bar + "background" + bar + "beam_data.mat",'ab') as f:
    #    scipy.io.savemat(f, {'K': K, 'I': I, 'M': M, 'C': C}) 
    
    # Simulation timestep:
    dt = model_param.timestep
    t0 = model_param.tspan[0]
    t1 = model_param.tspan[1]
    tvec = np.array([np.linspace(t0,t1,int((t1-t0)/dt)+1)])
    teval = matlab.double(tvec)
    
    # Initial condition parameters:
    np.random.seed(seed)
    p0span = model_param.p0span
    q0span = model_param.q0span
    
    # Dataset creation loop:
    X_data, Xdot_data, t_data = [], [], []
    for i in range(samples):
        p0 = (p0span[1]-p0span[0])*np.random.rand(n_DoF) + p0span[0] 
        q0 = (q0span[1]-q0span[0])*np.random.rand(n_DoF) + q0span[0]
        X0 = np.array([np.concatenate([p0,q0])])
        X0 = matlab.double(X0)
        
        X, Xdot = eng.simulation(model,DoF,gravity,teval,X0,disp_progress,nargout=2)
        
        X = np.array(X)        
        X_data.append(X)
        
        Xdot = np.array(Xdot)        
        Xdot_data.append(Xdot)
             
        t_data.append(tvec.T)
        
        progress_msg = '\rCreating dataset... {:.1f}% done'.format(100*(i+1)/samples)
        print(progress_msg + '\n' if i == samples-1 else progress_msg, end='')
    x_aux = np.concatenate(X_data)
    t = np.zeros_like(x_aux[...,:1])
    data_raw = {'x': np.concatenate(X_data), 
                'dx': np.concatenate(Xdot_data), 
                't': t}
    
    data_raw.update({'n_DoF': n_DoF, 'M': M, 'I': I, 'K': K, 'C': C})
    
    return data_raw
#data,data_raw = make_dataset(samples=2)