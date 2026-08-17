#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 09:13:20 2022

@author: dsw002
"""

# This is the same as lmfit_routines.py
# The only differences are 1) where the raw_data lives,
# and 2) a line I will designate which reframes the "raw"
# data as a random bootstrap sample.
# Call this function many times to get bootstrapped parameter estimates.

import os
import sys 
sys.path.append(os.getcwd())
sys.path.append('/home/gddaslab/dsw002/nk_cellular_development/clonal_expansion_python/for_publication/')
import numpy as np
import pandas as pd
import lmfit
import cost_functions as cf
from fit_moments_and_such import calc_moments_from_sample,bootstrap_sample
from clonal_burst_script import clonal_burst_2D as cb

def fit_moments_and_correlation(constraint):
    raw_data = pd.read_csv('../../../counts_CD27.csv')
    raw_data = raw_data.values
    ### HERE IS THE BOOTSTRAP CHANGE ###
    raw_data = bootstrap_sample(raw_data)
        
    bs_data = []
    bs_raw_data = []
    for i in range(10000):
        bs_sample = bootstrap_sample(raw_data)
        bs_data.append(calc_moments_from_sample(bs_sample))
        bs_raw_data.append(bs_sample)
    errors = np.std(bs_data,axis=0)
    data = calc_moments_from_sample(raw_data)
    
    params = lmfit.Parameters()
    params.add_many(('b1',1.06357,True,0,2),
                    ('b2',1.61259,True,0,2),
                    ('b3',0.69969, True,0,2),
                    ('d1',0.054277,True,0,2),
                    ('d2',0.28539,True,0,2),
                    ('d3',1.370,True,0,2),
                    ('diff1',0.1260,True,0,2),
                    ('diff2',0.16439,True,0,2))
    results = lmfit.minimize(cf.moment_residuals_3D_for_lmfit,params,args=(data,errors,True,False,constraint),method='least_squares')
    
    return results

if __name__ == '__main__':
    constraint = float(sys.argv[1])
    results = fit_moments_and_correlation(constraint)
    name = "output"+sys.argv[2]
    with open(name,"w") as fid:
        fid.write(str(results.params['b1'].value)+' '+str(results.params['b2'].value)+' '+str(results.params['b3'].value)+' '+str(results.params['d1'].value)+' '+str(results.params['d2'].value)+' '+str(results.params['d3'].value)+' '+str(results.params['diff1'].value)+' '    +str(results.params['diff2'].value)+' '+str(results.chisqr))
