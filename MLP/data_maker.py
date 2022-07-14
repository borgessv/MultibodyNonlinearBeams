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
import pandas as pd
import numpy as np
import matlab.engine

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

beam_data_file = parent_dir + "\\beam_data_test.ods"
beam_data = pd.read_excel(beam_data_file, engine="odf", skiprows=[0])

n_element = beam_data[1][2]

def make_dataset(seed=0, samples=20, test_split=0.75, **kwargs):
    print('\rCreating dataset... {:.1f}% done'.format(0), end='')
    eng = matlab.engine.start_matlab()
    #eng.create_beam("beam_data_test.ods",nargout=0)
    #data = {'meta': locals()}
    #np.random.seed(seed)
    dt = 0.1
    tf = 10
    tspan = [x*dt for x in range(int(tf/dt + 1))]
    
    q0max = np.pi
    p0max = 2*np.pi
    #q0_initial = q0span[0]*(2*np.random.rand(n_element)-1)
    #q0_final = q0_initial + (q0span[1]-q0span[0])*(2*np.random.rand()-1)
    X_data, Xdot_data, H_data = [], [], []

    for i in range(samples):
        
        q0 = q0max*(2*np.random.rand(n_element)-1) #np.linspace(q0_initial,q0_final,n_element,endpoint=True)[np.newaxis].T
        p0 = p0max*(2*np.random.rand(n_element)-1)
        X0 = np.concatenate([p0,q0])

        out = eng.simulation('FOM',{'OutBend'},'GravityOn',matlab.double(tspan),matlab.double(X0),nargout=3)

        p, q = np.split(np.array(out[0]).T,2)
        
        X_data.append(np.concatenate([q,p]).T)
        
        pdot, qdot = np.split(np.array(out[1]),2)
        Xdot_data.append(np.concatenate([qdot,pdot]).T)

        H = np.array(out[2])
        H_data.append(H)
        
        progress_msg = '\rCreating dataset... {:.1f}% done'.format(100*(i+1)/samples)
        print(progress_msg + '\n' if i == samples-1 else progress_msg, end='')
        
    
    data = {'X': np.concatenate(np.array(X_data)), 'dX': np.concatenate(np.array(Xdot_data))}
    #data['X'] = np.concatenate(X_data)
    #data['dX'] = np.concatenate(Xdot_data)

    # make a train/test split
    split_ix = int(len(data['X'])*test_split)
    split_data = {}
    for k in ['X', 'dX']:
        split_data[k], split_data['test_' + k] = data[k][:split_ix], data[k][split_ix:]
    data = split_data
    
    return data

#data = make_dataset(samples=2)