#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 10:40:03 2022

@author: dsw002
"""
# My Gillespie function.  Very robust to different inputs.

import numpy as np
import random

# inputs are:
# -init: the initial population (usually a CD27+ cell)
# -M: the kinetics matrix (not including death and asymmetric divison)
# -asym_M: Kinetics matrix describing how cells asymetrically divide
# -death: A vector of death rates for each cell type
# -end_t: The time we count cells at

# An example input:
# For a 3-stage system with these parameters:
# birth rates: 0.3, 0.2, 0.1
# death rates: 0.1, 0.1, 0
# differentiation rates: 0.5, 0.4
# asymmetric division rate of 1st to 2nd stage: 0.1
# init = [1,0,0]
# M = [[0.3, 0, 0],[0.5, 0.2, 0],[0, 0.4, 0.1]]
# asym_M = [[0.1, 0, 0],[0.1, 0, 0],[0, 0, 0]]
# death = [0.1, 0.1, 0]
# end_t = 8

def gillespie_for_asymmetric_division(init,M,asym_M,death,end_t):
    # initialize the algorithm
    n_species = init
    t = 0
    continue_bool = True
    # until we reach the end_t
    while continue_bool:
        # determine propensities for each reaction given reaction rates and abundances
        propensities = np.row_stack((np.tile(n_species,[2*len(n_species),1]), n_species)) * np.abs(np.row_stack((M, asym_M, death)))
        # determine the time step that will occur
        delta_t = np.log(1/random.random())*(1/np.sum(propensities))
        # transform propensities so we can determine which reaction to choose
        propensities = np.cumsum(propensities/np.sum(propensities))
        t = t+delta_t
        # if we haven't reached the end
        if t<end_t:
            # select a reaction at random
            rxn_rand = random.random()
            which_rxn = next(i for i,val in enumerate(propensities) if rxn_rand<val)+1
            # if that reaction is in the first n^2, that means it is not asym division
            # or death
            if which_rxn<=len(n_species)**2:
                # the "sink" and "source" describe the product and reactant, respectively
                sink = int(np.ceil(which_rxn/len(n_species)))-1
                source = int(np.mod(which_rxn,len(n_species)))-1
                # Due to the modulo operations, if the reactant is the last species
                # we need to reset it to be the index of the last species
                if source==-1:
                    source = len(n_species)-1
                # if the reaction is on the diagonal, meaning it is a proliferation event
                if source==sink:
                    n_species[source] = n_species[source]+1
                else: # otherwise, it's differentiation (conversion)
                    n_species[sink] = n_species[sink]+1
                    n_species[source] = n_species[source]-1
            
            # if the reaction is a asymmetric division one
            elif which_rxn<=2*len(n_species)**2:
                # update the value so we can repeat the process we did in the previous 'if'
                # statement
                which_rxn = which_rxn-len(n_species)**2
                sink = int(np.ceil(which_rxn/len(n_species)))-1
                source = int(np.mod(which_rxn,len(n_species)))-1
                if source==-1:
                    source = len(n_species)-1
                if source==sink: #in theory, this should never occur
                    n_species[source] = n_species[source]+1
                else:
                    # instead, generate a cell as a product
                    n_species[sink] = n_species[sink]+1
            
            # if not proliferation, differentiation, or asym division, it's death
            else:
                n_species[which_rxn-2*len(n_species)**2-1] = n_species[which_rxn-2*len(n_species)**2-1]-1
        else: # if we reached the end_t
            continue_bool = False
            
    return n_species

if __name__ == "__main__":
    cells = gillespie_for_asymmetric_division(np.array([1,0]),np.array([[1.1, 0],[0.5,0.2]]),np.array([[0,0],[0,0]]),np.array([0,0]),8)
