#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 09:13:20 2022

@author: dsw002
"""

# This code is a "wrapper" to fit the moments of the NK data given a set of parameters.
# Here is where you will define your system, initial guess, and constraints.

import os
import sys 
sys.path.append(os.getcwd())
sys.path.append('/home/gddaslab/dsw002/nk_cellular_development/clonal_expansion_python/for_publication/')
import numpy as np
import pandas as pd
import lmfit

# These are modules I wrote:
import cost_functions as cf
from fit_moments_and_such import calc_moments_from_sample,bootstrap_sample
from clonal_burst_script import clonal_burst_2D as cb

# My fitting routine
def fit_moments_and_correlation(constraint):
    # Read the data
    raw_data = pd.read_csv('counts_CD27.csv')
    raw_data = raw_data.values
    
    # Get bootstrapped data to get uncertainty in moments
    np.random.seed(12345)
    bs_data = []
    bs_raw_data = []
    for i in range(10000): #10000 bootstrap samples
        # get a sample "experiment"
        bs_sample = bootstrap_sample(raw_data)
        # calculate its moments
        bs_data.append(calc_moments_from_sample(bs_sample))
        # hold onto that data in case we need it (it doesn't look like I ever used this)
        bs_raw_data.append(bs_sample)
    # uncertainty in my data's moments are standard deviations of the bootstrapped moments
    errors = np.std(bs_data,axis=0)
    # get my data's moments
    data = calc_moments_from_sample(raw_data)
    print(data)
    print(errors)
    # reset seed in case stochastic optimizer is used
    np.random.seed()
    
    # set initial parameter estimates. This is the area to change for varying systems.
    params = lmfit.Parameters()
    # Each entry is (parameter_name,initial_value,fit_it?,minimum,maximum)
    # I HIGHLY recommend checking out the lmfit documentation.  It's very good and will be useful for the future.
    params.add_many(('b1',0.64203,True,0,2),
                    ('b2',1.42808,True,0,2),
                    ('b3',0.699087, True,0,2),
                    ('d1',0.00755127,True,0,2),
                    ('d2',0.422314,True,0,2),
                    ('d3',1.38265,True,0,2),
                    ('diff1',0.248496,True,0,2),
                    ('diff2',0.260533,True,0,2))
    # call the fitter.  Arguments are: (moments,errors/uncertainty,is_3_stage?,is_2_mature?,correlation_constraint)
    results = lmfit.minimize(cf.moment_residuals_3D_for_lmfit,params,args=(data,errors,True,False,constraint),method='least_squares')
    
    return results

if __name__ == '__main__':
    # The first argument when calling the script is the constraint for the correlation (i.e. -0.2)
    constraint = float(sys.argv[1])
    results = fit_moments_and_correlation(constraint)
    # write output
    with open("a0_result","w") as fid:
        fid.write(str(results.params['b1'].value)+' '+str(results.params['b2'].value)+' '+str(results.params['b3'].value)+' '+str(results.params['d1'].value)+' '+str(results.params['d2'].value)+' '+str(results.params['d3'].value)+' '+str(results.params['diff1'].value)+' '    +str(results.params['diff2'].value)+' '+str(results.chisqr))
