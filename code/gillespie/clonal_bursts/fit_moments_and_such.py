#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 09:29:18 2022

@author: dsw002
"""
# Much of this is depracated, but some functions still hold value to other scripts.
# I've denoted where these are.  Notably, calculation of moments and bootstrapping
# samples are done here.

import numpy as np
from scipy.linalg import expm
import lmfit
import matplotlib.pyplot as plt
import pandas as pd

# calculate the moments of a dataset
def calc_moments_from_sample(data):
    n = np.mean(data[:,0])
    m = np.mean(data[:,1])
    nn = np.mean(data[:,0]*data[:,0])
    mm = np.mean(data[:,1]*data[:,1])
    nm = np.mean(data[:,0]*data[:,1])
    return [n,m,nn,mm,nm]

# deprecated; see cost_functions.py for updated
def moment_solutions(b1,b2,d1,d2,diff,asym,t):
    M = np.array([[b1-d1-diff,0,0,0,0],[diff+asym,b2-d2,0,0,0],[b1+d1+diff,0,2*(b1-d1-diff),0,0],[diff+asym,b2+d2,0,2*(b2-d2),2*(diff+asym)],[-diff,0,diff+asym,0,b1-d1+b2-d2-diff]])
    init = [1,0,1,0,0]
    solutions = np.dot(expm(M*t),init)
    return solutions

# deprecated; see cost_functions.py for updated
def calc_residual(params,data,errors):
    fit = moment_solutions(params['b1'],params['b2'],params['d1'],params['d2'],params['diff'],params['asym'],8)
    resid = (data-fit)/errors
    if params['b1'].value-params['d1'].value<params['b2'].value-params['d2'].value:
        boundary = 1000.
    else:
        boundary = 0.
    return np.append(resid,boundary)
    #cyt_resid = cytof_residuals(params)
    #return np.concatenate([resid,cyt_resid])

# deprecated; see cost_functions.py for updated
def cytof_residuals(params):
    data = np.array([[ 707677.34621169, 1273897.6167672 ],
    [ 319545.75551097, 4599077.83595458]])
    std_data = np.array([[46187.00009816, 95949.13237813],[42718.4756798, 46990.0702212]])
    M = np.array([[params['b1']-params['d1']-params['diff']+params['asym'],0],[params['diff']+params['asym'],params['b2']-params['d2']]])
    times = [4.0,7.0]
    resid = np.zeros((2,len(times)))
    for ind,t in enumerate(times):
        resid[:,ind] = (data[ind,:] - np.dot(expm(M*t),[307581,666482]))/std_data[ind,:]
    return resid.flatten()

# Return a random sample of the data for bootstrapping
def bootstrap_sample(data):
    # pick random indices of data
    ind = np.random.choice(range(len(data)),len(data))
    # return data according to the samples
    bs_sample = data[ind,:]
    return bs_sample

# Given a set of parameters, see how well it fits CyTOF data.
# I don't think I use this in the pipeline, and the data should be validated in one of the
# other python files before use.
def compare_to_cytof(params):
    M = np.array([[params['b1']-params['d1']-params['diff'],0],[params['diff'],params['b2']-params['d2']]])
    times = np.linspace(0,7,50)
    solution = np.zeros((2,len(times)))
    for ind,t in enumerate(times):
        solution[:,ind] = np.dot(expm(M*t),[307581,666482])
    plt.figure()
    plt.semilogy([0,4,7],np.array([[ 307581.26107153,  666482.10893078],
           [ 707677.34621169, 1273897.6167672 ],
           [ 319545.75551097, 4599077.83595458]]),'o')
    plt.gca().set_prop_cycle(None)
    plt.semilogy(times,solution.T)
    plt.xlabel('Days')
    plt.ylabel('#Cells')
    plt.legend(['CD27+','CD27-'])
    
# If passed bootstrapped parameter values, get confidence intervals.
# a0_case is the parameter set that best fits the data, params is a list
# of bootstrapped parameter sets.
# The output is a 3xn matrix, with the lower_bound, estimate, and upper_bound
# in each row respectively.
# I think I have an updated version of this, as this only works for a 
# 2-stage model with asymmetric division.  Probably deprecated.
def get_outputs(params,a0_case):
    
    # process the input
    d1 = []
    d2 = []
    k1 = []
    k2 = []
    k3 = []
    asym = []
    for p in params:
        d1.append(p['b1'])
        d2.append(p['b2'])
        k1.append(p['d1'])
        k2.append(p['d2'])
        k3.append(p['diff'])
        asym.append(p['asym'])
    
    # determine confidence intervals by figuring out 2.5th and 97.5th percentile
    # values from bootstrapped params
    outputs = []
    outputs.append((np.percentile(d1,2.5),a0_case.params['b1'].value,np.percentile(d1,97.5)))
    outputs.append((np.percentile(d2,2.5),a0_case.params['b2'].value,np.percentile(d2,97.5)))
    outputs.append((np.percentile(k1,2.5),a0_case.params['d1'].value,np.percentile(k1,97.5)))
    outputs.append((np.percentile(k2,2.5),a0_case.params['d2'].value,np.percentile(k2,97.5)))
    outputs.append((np.percentile(k3,2.5),a0_case.params['diff'].value,np.percentile(k3,97.5)))
    outputs.append((np.percentile(asym,2.5),a0_case.params['asym'].value,np.percentile(asym,97.5)))
    return np.asarray(outputs)

# deprecated; see lmfit_routines.py for updated
def main():
    raw_data = pd.read_csv('~/Dropbox/Work/NK Cell Development/clonal_expansion_python/counts_CD27.csv')
    raw_data = raw_data.values
    
    bs_data = []
    bs_raw_data = []
    for i in range(10000):
        bs_sample = bootstrap_sample(raw_data)
        bs_data.append(calc_moments_from_sample(bs_sample))
        bs_raw_data.append(bs_sample)
    errors = np.std(bs_data,axis=0)

    data = calc_moments_from_sample(raw_data)
    params = lmfit.Parameters()
    params.add_many(('b1',2,True,0.1,10),
                    ('b2',1,True,0.1,10),
                    ('d1',1,True,0.01,10),
                    ('d2',1,True,0.01,10),
                    ('diff',1,True,0.001,10),
                    ('asym',0,False))
    a0_results = lmfit.minimize(calc_residual,params,args=(data,errors))

    plt.figure()
    plt.scatter(np.log10(data),np.log10(data-a0_results.residual[0:5]*errors))
    plt.plot(np.log10([np.min(data),np.max(data)]),np.log10([np.min(data),np.max(data)]))
    
    compare_to_cytof(a0_results.params)

    bs_params = []
    for bs_sample in bs_raw_data:
        second_bs_data = []
        for i in range(10000):
            second_bs_sample = bootstrap_sample(bs_sample)
            second_bs_data.append(calc_moments_from_sample(second_bs_sample))
        bs_errors = np.std(second_bs_data,axis=0)
        in_data = calc_moments_from_sample(bs_sample)
        result = lmfit.minimize(calc_residual,params,args=(in_data,bs_errors))
        bs_params.append(result.params)
    
    cis = get_outputs(bs_params,a0_results)
    return cis
        

if __name__ == "__main__":
    cis = main()



