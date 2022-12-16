# -*- coding: utf-8 -*-
"""
"""
import torch
import scipy.integrate
import numpy as np
import os, sys
import matlab.engine

if 'windows' in sys.platform:
    bar = '\\'
else:
    bar = '/'

def matlab_interface():
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(parent_dir) 
    eng = matlab.engine.start_matlab()
    eng.addpath(parent_dir + bar + "SimulationFramework", nargout= 0)
    eng.addpath(parent_dir + bar + "SimulationFramework" + bar + "background", nargout= 0)
    eng.addpath(parent_dir + bar + "SimulationFramework" + bar + "background" + bar + "utils", nargout= 0)
    eng.addpath(parent_dir + bar + "SimulationFramework" + bar + "background" + bar + "CrossSectionData", nargout= 0)
    return eng


# simplifies accessing the hyperparameters.
class ObjectView(object):
    def __init__(self, d): self.__dict__ = d
    

def integrate_model(model, t_span, y0, n_DoF, use_torch=True, **kwargs):
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


def normalize(data, norm_range):
    n_DoF = int(data['x'].shape[1]/2)
    data['x'][:,:n_DoF] = np.ptp(norm_range)*(data['x'][:,:n_DoF] - np.min(data['x'][:,:n_DoF]))/np.ptp(data['x'][:,:n_DoF]) + norm_range[0]   
    data['x'][:,n_DoF:] = np.ptp(norm_range)*(data['x'][:,n_DoF:] - np.min(data['x'][:,n_DoF:]))/np.ptp(data['x'][:,n_DoF:]) + norm_range[0]
    data['dx'][:,:n_DoF] = np.ptp(norm_range)*(data['dx'][:,:n_DoF] - np.min(data['dx'][:,:n_DoF]))/np.ptp(data['dx'][:,:n_DoF]) + norm_range[0]    
    data['dx'][:,n_DoF:] = np.ptp(norm_range)*(data['dx'][:,n_DoF:] - np.min(data['dx'][:,n_DoF:]))/np.ptp(data['dx'][:,n_DoF:]) + norm_range[0]
    return data


def standardize(data):
    n_DoF = int(data['x'].shape[1]/2)
    data_new = data
    data_new['x'][:,:n_DoF] = (data['x'][:,:n_DoF] - np.mean(data['x'][:,:n_DoF]))/np.std(data['x'][:,:n_DoF])   
    data_new['x'][:,n_DoF:] = (data['x'][:,n_DoF:] - np.mean(data['x'][:,n_DoF:]))/np.std(data['x'][:,n_DoF:])
    data_new['dx'][:,:n_DoF] = (data['dx'][:,:n_DoF] - np.mean(data['dx'][:,:n_DoF]))/np.std(data['dx'][:,:n_DoF])    
    data_new['dx'][:,n_DoF:] = (data['dx'][:,n_DoF:] - np.mean(data['dx'][:,n_DoF:]))/np.std(data['dx'][:,n_DoF:])
    return data_new
    
    