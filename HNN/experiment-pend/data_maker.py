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

def make_dataset(seed=0, samples=20, test_split=0.5, **kwargs):
    print('\rCreating dataset... {:.1f}% done'.format(0), end='')
    
    # Model parameters:    
    model = "FOM"
    DoF = {"OutBend"}
    gravity = "GravityOn"
    
    # Initializing structure properties:
    beam_data = "beam_data_test.xlsx"
    beam = pd.read_excel(parent_dir + "\\SimulationFramework\\" + beam_data)
    n_DoF = len(DoF)*int(np.sum(np.array(beam.iloc[3,1:])))
    M, I, K, C = eng.structure_properties(beam_data,DoF,nargout=4)
    # K = eng.times(1e-2,K)
    K = matlab.double(K)
    K[0][0] = 0
    C = eng.times(5e-5,C)
    
    # Simulation timestep:
    dt = 0.1
    tf = 10
    tspan = [x*dt for x in range(int(tf/dt + 1))]
    tspan = matlab.double(tspan)
    
    # Initial condition parameters:
    np.random.seed(seed)
    q0max = np.pi
    p0max = 0
    
    # Dataset creation loop:
    X_data, Xdot_data = [], []
    for i in range(samples):
        q0 = q0max*(2*np.random.rand(n_DoF)-1) 
        p0 = p0max*(2*np.random.rand(n_DoF)-1)
        X0 = np.array([np.concatenate([p0,q0])]).T
        X0 = matlab.double(X0)
        
        # #if y0 is None:
        # y0 = np.random.rand(2)*2.-1
        # #if radius is None:
        # radius = np.random.rand() + 1.3 # sample a range of radii
        # y0 = y0/np.sqrt((y0**2).sum())*radius ## set the appropriate radius
        # X0 = matlab.double(y0)
        
        X, Xdot = eng.simulation(model,DoF,gravity,tspan,M,I,K,C,X0,nargout=2)
        
        X = np.array(X)        
        X_data.append(X)
        
        Xdot = np.array(Xdot)        
        Xdot_data.append(Xdot)
        
        progress_msg = '\rCreating dataset... {:.1f}% done'.format(100*(i+1)/samples)
        print(progress_msg + '\n' if i == samples-1 else progress_msg, end='')
           
    data = {'X': np.concatenate(np.array(X_data)), 'dX': np.concatenate(np.array(Xdot_data))}

    # make a train/test split
    split_ix = int(len(data['X'])*test_split)
    split_data = {}
    for k in ['X', 'dX']:
        split_data[k], split_data['test_' + k] = data[k][:split_ix], data[k][split_ix:]
    data = split_data
    
    return data

# data = make_dataset(samples=2)