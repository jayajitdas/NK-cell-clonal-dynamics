#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 09:36:08 2022

@author: dsw002
"""
# gillespie script is my own
from gillespie_for_asymmetric_division import gillespie_for_asymmetric_division
import numpy as np
import numpy.ma as ma

# Note: the correlation is calculated with masked_invalid, this is because NAN values are possible and will ruin the calculation without it

# Clonal burst experimental correlation with 2 cell types
def clonal_burst_2D(b1,b2,d1,d2,diff,asym,n_iters):
    # initialize the clones
    cells_at_t = np.zeros((n_iters,2))
    for i in range(n_iters): #simulate each clone
        cells_at_t[i,:]=gillespie_for_asymmetric_division(np.array([1,0]),np.array([[b1, 0],[diff,b2]]),np.array([[0,0],[asym,0]]),np.array([d1,d2]),8)
    # calculate burst size and %CD27 for each clone
    burst_size = np.sum(cells_at_t,axis=1)
    cd27_percent = cells_at_t[:,0]/burst_size
    # calculate correlation between those values
    corr = ma.corrcoef(ma.masked_invalid(cd27_percent), ma.masked_invalid(burst_size))
    return corr[0,1]

# Clonal burst experimental correlation with 3 cell types
def clonal_burst_3D(b1,b2,b3,d1,d2,d3,diff1,diff2,n_iters,one_mature=True):
    cells_at_t = np.zeros((n_iters,3))
    for i in range(n_iters):
        # How should initial condition be determined? If "one_mature", 73% of CD27+ cells are in the second stage, and 27% are in the first
        if one_mature:
            if np.random.rand()<0.2713:
                init = np.array([1,0,0])
            else:
                init = np.array([0,1,0])
        # If two mature stages, just start with the first stage
        else:
            init = np.array([1,0,0])
        # simulate each clone
        cells_at_t[i,:]=gillespie_for_asymmetric_division(init,np.array([[b1, 0, 0],[diff1,b2,0],[0,diff2,b3]]),np.array([[0,0,0],[0,0,0],[0,0,0]]),np.array([d1,d2,d3]),8)
    # calculate burst size and %CD27
    burst_size = np.sum(cells_at_t,axis=1)
    if one_mature:
        cd27_percent = (cells_at_t[:,0]+cells_at_t[:,1])/burst_size
    else:
        cd27_percent = (cells_at_t[:,0])/burst_size
    # calculate correlation between these
    corr = ma.corrcoef(ma.masked_invalid(cd27_percent), ma.masked_invalid(burst_size))
    return corr[0,1]

# This shouldn't normally be used; this is more for test or validation
if __name__ == "__main__":
    corr = clonal_burst_3D(0.2,0.8,0.25,0,0.2,0.2,0.1,0.5,1000)
    print(corr)
