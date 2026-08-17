#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 11:37:14 2022

@author: dsw002
"""
# This is a set of functions that can be used to fit cost functions to data
# NOTE: Here is where you will want to make edits for new model fitting routines
# I'll notate them, but see "determine_pdf" in particular.


import numpy as np
from scipy.linalg import expm #matrix exponentiation
# I wrote the following modules
from clonal_burst_script import clonal_burst_2D as cb
from clonal_burst_script import clonal_burst_3D as cb3

# Calculate likelihood given a set of parameters (that will generate a PDF) and data
def KDE_likelihood(params,data):
    # generate the pdf
    model_pdf = determine_pdf(params)
    # Calculate the probability. This is a product of all probabilities, so start with p=1
    p = 1
    # multiply each probability together
    for point in data:
        # the probability of each point is given by the value of the pdf at the data point
        p = p*model_pdf(data[0],data[1])
    # cost is -2 log likelihood
    cost = -2*np.log(p)
    return cost

#### USE THIS FOR GENERATING PDF ####
# The output should be a matrix with each value as a discrete value of number of immature
# or mature cells.
def determine_pdf(params):
    
    return pdf

# Calculate residuals from CyTOF data if only gating CD27+ vs. CD27-
def cytof_residuals_2D(b1,b2,d1,d2,diff,asym):
    # hardcoded data and uncertainty
    data = np.array([[ 707677.34621169, 1273897.6167672 ],
    [ 319545.75551097, 4599077.83595458]])
    std_data = np.array([[46187.00009816, 95949.13237813],[42718.4756798, 46990.0702212]])
    # kinetics matrix
    M = np.array([[b1-d1-diff+asym,0],[diff+asym,b2-d2]])
    # times post-infection CyTOF was measured at
    times = [4.0,7.0]
    # calculate residuals for each time
    resid = np.zeros((2,len(times)))
    for ind,t in enumerate(times):
        resid[:,ind] = (data[ind,:] - np.dot(expm(M*t),[307581,666482]))/std_data[ind,:]
    return resid

# moment solutions for "2D" system, meaning only CD27+ and CD27- cells are gated
def moment_solutions(b1,b2,d1,d2,diff,asym,t):
    # moments are [immature, mature, immature^2, mature^2, immature*mature]
    
    # first order kinetics of moment solution
    M = np.array([[b1-d1-diff,0,0,0,0],[diff+asym,b2-d2,0,0,0],[b1+d1+diff,0,2*(b1-d1-diff),0,0],[diff+asym,b2+d2,0,2*(b2-d2),2*(diff+asym)],[-diff,0,diff+asym,0,b1-d1+b2-d2-diff]])
    init = [1,0,1,0,0]
    solutions = np.dot(expm(M*t),init)
    return solutions

# residuals given calculated moments (from parameters) and data
def moment_residuals(b1,b2,d1,d2,diff,asym,data,errors):
    fit = moment_solutions(b1,b2,d1,d2,diff,asym,8)
    resid = (data-fit)/errors
    return resid

# updated moment_residuals code such that it can handle the Parameters() input.
# It also now has additional constraints that can be passed to it as input args.
def moment_residuals_for_lmfit(params,data,errors,corr_bound=True,constraint=1):
    fit = moment_solutions(params['b1'],params['b2'],params['d1'],params['d2'],params['diff'],params['asym'],8)
    resid = (data-fit)/errors
    # If immature growth is slower, add an arbitrary cost
    if (params['b1']-params['d1']+params['asym'])<(params['b2']-params['d2']):
        resid = np.concatenate((resid,[10000]))
    else: # if not, add nothing
        resid = np.concatenate((resid,[0]))
    # if a correlation constraint is provided
    if corr_bound:
        # evaluate the correlation from Gillespie
        corr = cb(params['b1'],params['b2'],params['d1'],params['d2'],params['diff'],params['asym'],1000)
        print(corr)
        # if the correlation is not sufficiently negative, add an arbitrary cost
        if corr>constraint:
            resid = np.concatenate((resid,[10000]))
        else:
            resid = np.concatenate((resid,[0]))
    return resid

# Same as moment_solutions, but for data that is gated into 3 populations rather than 2
# input argument "one_mature" designates whether the intermediate stage is mature or not.
# If "True", then it is immature, if "False" then it is mature.
def moment_solutions_3D(b1,b2,b3,d1,d2,d3,diff1,diff2,t=8,one_mature=True):
    #3 cell types: n->m->l
    #moments: [n,m,l,n^2,m^2,l^2,nm,ml,nl]
    M = np.array([[b1-d1-diff1,0,0,0,0,0,0,0,0],
                [diff1,b2-d2-diff2,0,0,0,0,0,0,0],
                [0,diff2,b3-d3,0,0,0,0,0,0],
                [b1+d1+diff1,0,0,2*(b1-d1-diff1),0,0,0,0,0],
                [diff1,b2+d2+diff2,0,0,2*(b2-d2-diff2),0,2*diff1,0,0],
                [0,diff2,b3+d3,0,0,2*(b3-d3),0,2*diff2,0],
                [-diff1,0,0,diff1,0,0,b1+b2-d1-d2-diff1-diff2,0,0],
                [0,-diff2,0,0,diff2,0,0,b2+b3-d2-d3-diff2,diff1],
                [0,0,0,0,0,0,diff2,0,b1+b3-d1-d3-diff1]])
    # This system needs special initial conditions.
    # I've chosen the case where the 2nd stage is immature to reflect the CD11b/CD27
    # CyTOF gating proportions.
    # If the 2nd stage is mature, I'll initialize only with the 1st stage.
    if one_mature:
        init = [0.2713,0.7287,0,0.2713,0.7287,0,0,0,0]
    else:
        init = [1,0,0,1,0,0,0,0,0]
    
    solutions = np.dot(expm(M*t),init)
    # calculate the moments for CD27+ and CD27- cells as a whole
    if one_mature:
        EI = solutions[0]+solutions[1] #E[n]+E[m]
        EM = solutions[2] #E[l]
        EI2 = solutions[3]+solutions[4]+2*solutions[6] #E[n^2]+E[m^2]+2E[nm]
        EM2 = solutions[5] #E[l^2]
        EIM = solutions[7]+solutions[8] #E[ml]+E[nl]
    else:
        EI = solutions[0]
        EM = solutions[1]+solutions[2]
        EI2 = solutions[3]
        EM2 = solutions[4]+solutions[5]+2*solutions[7]
        EIM = solutions[6]+solutions[8]
    return np.array([EI,EM,EI2,EM2,EIM])

# Same as moment_residuals except for 3 NK maturation stages
def moment_residuals_3D(b1,b2,b3,d1,d2,d3,diff1,diff2,data,errors,corr_bound=True,one_mature=True):
    fit = moment_solutions_3D(b1,b2,b3,d1,d2,d3,diff1,diff2,8,one_mature)
    resid = (data-fit)/errors
    return resid

# Same as moment_residuals_for_lmfit except for 3 NK maturation stages
def moment_residuals_3D_for_lmfit(params,data,errors,corr_bound=True,one_mature=True,constraint=1):
    fit = moment_solutions_3D(params['b1'],params['b2'],params['b3'],params['d1'],params['d2'],params['d3'],params['diff1'],params['diff2'],8,one_mature)
    resid = (data-fit)/errors
    M = np.array([[params['b1']-params['d1'],0,0],[params['diff1'],params['b2']-params['d2'],0],[0,params['diff2'],params['b3']-params['d3']]])
    # unlike for the 2-stage case, we can't compare parameters to determine if CD27+
    # cells grow faster than CD27-.  Thus we need to do a calculation for it. This does
    # that, with different "effective" growth rates depending on the initial condition.
    # subset_1_bool indicates whether the cell type is CD27+.
    if not(one_mature):
        init = np.array([1,1,1])
        subset_1_bool = np.array([True,False,False])
    else:
        init = np.array([0.08223464,0.2208308275,1])
        subset_1_bool = np.array([True,True,False])
    # calculate the difference in observed NK cell populations resulting from CD27+ and 
    # CD27- cells.
    difference = difference_for_CD27pos(M,init,subset_1_bool,False)
    # if CD27- cells grow faster, add an arbitrary cost
    if difference<0:
        resid = np.concatenate((resid,[10000]))
    else:
        resid = np.concatenate((resid,[0]))
    if corr_bound:
        corr = cb3(params['b1'],params['b2'],params['b3'],params['d1'],params['d2'],params['d3'],params['diff1'],params['diff2'],10000,one_mature)
        print(corr)
        if corr>constraint:
            resid = np.concatenate((resid,[10000]))
        else:
            resid = np.concatenate((resid,[0]))
    return resid
  
# Calculates the number of cells that would result from an initial population of CD27+ 
# cells, and compares that to the number that would result from an initial population
# of CD27- cells, and returns the difference (#CD27+ - #CD27-)
# 'subtracts_diff' is some functionality that I can specify whether I gave it just the
# rates as a matrix or whether I gave the kinetics matrix (i.e. are the diagonals the 
# growth rates, or the growth rates minus the sum of the column? If they are just the
# rates, then subtracts_diff is False)
def difference_for_CD27pos(M,initial_abundance,subset_1_bool,subtracts_diff=False):
    if not(subtracts_diff):
        for i in range(len(M[0,:])):
            M[i,i] = 2*M[i,i]-np.sum(M[:,i])
    n_cells = np.zeros((len(M[0,:]),))
    # for each cell type, how many cells would result from a single cell of that type
    for i in range(len(M[0,:])):
        init = np.zeros((len(M[0,:]),))
        init[i] = 1 # at this point, the initial condition should have one non-zero element
        # evaluate the kinetics
        n_cells[i] = np.sum(np.dot(expm(M*8),init))
    # determine what would result from the "average" CD27+ or CD27- cell
    initial_abundance[subset_1_bool] = initial_abundance[subset_1_bool]/np.sum(initial_abundance[subset_1_bool])
    initial_abundance[~subset_1_bool] = initial_abundance[~subset_1_bool]/np.sum(initial_abundance[~subset_1_bool])
    n_subset_1 = np.sum(np.log(n_cells[subset_1_bool])*initial_abundance[subset_1_bool])
    n_subset_2 = np.sum(np.log(n_cells[~subset_1_bool])*initial_abundance[~subset_1_bool])
    
    difference = n_subset_1-n_subset_2
    return difference
 
