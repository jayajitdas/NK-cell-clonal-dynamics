#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 15:05:26 2021##

@author: darrenwethington
"""

import numpy as np
import pandas as pd
import lmfit
from scipy.linalg import expm
import glob
import matplotlib.pyplot as plt

Ly49H_pos = 1

def model_process(params,time,x0,data,error):
    M = np.asarray([[params['k1'],0,0],
                    [params['d1'],params['k1'],0],
                    [0,params['d2'],params['k3']]])
    for i in range(np.size(M,axis=0)):
        M[i,i] = 2*M[i,i]-np.sum(M[:,i])
    y = np.zeros((np.size(time),np.size(data,axis=0)))
    for idx,t in enumerate(time):
        y[idx,:] = np.dot(expm(M*t),x0)
    resid = (np.log(y.T)-np.log(data))/error
    if np.any(np.isnan(resid)):
        print(resid)
        print(M)
    return resid
    

def call_model(data,time,error):
    data = np.asarray(data)
    params = lmfit.Parameters()
    params.add_many(('d1',0.1,True,0),
                    ('d2',0.1,True,0),
                    ('k1',0.1,True,0),
                    ('k2',0.1,True),
                    ('k3',0.1,True))
    #Order: [H,11b,27]. Ex: data[2] == H+,11b+,27
    results = lmfit.minimize(model_process,params,args=(time[1:],data[0,[2,3,1]].T,data[1:,[2,3,1]].T,error[1:,[2,3,1]].T))
    return results

def get_total_cell_counts():
    data = np.array([1014900.40000000,	2184353.38571429,	5017909,	961903.500000000,	1014900.40000000])
    return data

def read_data(file_string,idx,threshold):
    data = []
    percent = []
    #directory = '/Users/darrenwethington/Dropbox/Work/NK Cell Development/Work for Lanier/Data/new_data/MP311 Group 1 ILCs Only/Wild-Type Mice/'
    directory = '/Users/dsw002/Dropbox (NCH)/Work/NK Cell Development/Work for Lanier/Data/new_data/MP311 Group 1 ILCs Only/Wild-Type Mice/'
    for file in glob.glob(directory+'*.csv'):
        if file_string in file:
            x = pd.read_csv(file,index_col=None, header=0)
            #data.append(x)
            #percent_vec = np.zeros(np.size(idx)**2)
            boo = [False,True]
            full_boo = np.zeros(np.size(idx))
            percent_vec = []
            while full_boo[-1]<2:
                temp = np.ones(np.size(x,axis=0))
                for k,c in enumerate(idx):
                    temp = temp&((x.iloc[:,idx[k]]>threshold[k])==full_boo[k])
                
                full_boo[0] = full_boo[0]+1
                for i in range(np.size(full_boo)-1):
                    if full_boo[i]==2:
                        full_boo[i] = 0
                        full_boo[i+1] = full_boo[i+1]+1
                    
                percent_vec.append(np.average(temp))
            percent.append(percent_vec)
    return percent

def bootstrap_values(percent,Ly49H_pos):
    n=3
    #percent = np.asarray(percent)
    n_bootstrap = 10000
    for i in range(n):
        if Ly49H_pos:
            percent[i] = np.asarray(percent[i])[:,[1,3,5,7]]
        else:
            percent[i] = np.asarray(percent[i])[:,[0,2,4,6]]
        percent[i] = (percent[i].T/np.sum(percent[i],axis=1)).T
    percent_bootstrapped = np.zeros([n_bootstrap,4,n])
    percent_std = np.zeros([n_bootstrap,4,n])
    percent_pos = np.zeros((n,4))
    error_pos = np.zeros((n,4))
    for i in range(n):
        percent_pos[i,:] = np.average(percent[i],axis=0)
        error_pos[i,:] = np.std(percent[i],axis=0)
        for j in range(n_bootstrap):
            which_sample = np.random.choice((np.size(percent[i],axis=0)),(np.size(percent[i],axis=0)))
            bootstrap_sample = [percent[i][k] for k in which_sample]
            percent_bootstrapped[j,:,i] = np.average(bootstrap_sample,axis=0)
            percent_std[j,:,i] = np.std(bootstrap_sample,axis=0)
    return percent_bootstrapped,percent_pos,percent_std,error_pos

def get_outputs(params,a0_case):
    
    d1 = []
    d2 = []
    k1 = []
    k2 = []
    k3 = []
    for p in params:
        d1.append(p['d1'])
        d2.append(p['d2'])
        k1.append(p['k1'])
        k2.append(p['k2'])
        k3.append(p['k3'])
    
    
    outputs = []
    outputs.append((np.percentile(d1,2.5),a0_case.params['d1'].value,np.percentile(d1,97.5)))
    outputs.append((np.percentile(d2,2.5),a0_case.params['d2'].value,np.percentile(d2,97.5)))
    outputs.append((np.percentile(k1,2.5),a0_case.params['k1'].value,np.percentile(k1,97.5)))
    outputs.append((np.percentile(k2,2.5),a0_case.params['k2'].value,np.percentile(k2,97.5)))
    outputs.append((np.percentile(k3,2.5),a0_case.params['k3'].value,np.percentile(k3,97.5)))
    return np.asarray(outputs)

def get_outputs_that_satisfy_condition(params):
    d1 = []
    d2 = []
    k1 = []
    k2 = []
    k3 = []
    outputs = []
    for p in params:
        d1.append(p['d1'].value)
        d2.append(p['d2'].value)
        k1.append(p['k1'].value)
        k2.append(p['k2'].value)
        k3.append(p['k3'].value)
        if (p['k1'].value<p['k2'].value) & (p['k3'].value<p['k1'].value):
            outputs.append([p['k1'].value,p['k2'].value,p['k3'].value,p['d1'].value,p['d2'].value])
    #subset = (k1<k2) and (k3<k1)
    #print(subset)
    #outputs = [k1[subset],k2[subset],k3[subset],d1[subset],d2[subset]]
    return outputs


n = get_total_cell_counts()
n=n[0:3]
time = [0.0,4.0,7.0]
data = []
for idx,day in enumerate([0,4,7]):
    file_string = 'D'+str(day)
    data2 = read_data(file_string,[22,2,53],[np.sinh(3),np.sinh(2.5),0.0])
    data.append(data2)
percent_bootstrapped,percent_a0,error,error_a0 = bootstrap_values(data,Ly49H_pos)
a0_result = call_model((percent_a0.T*n).T,time,error_a0)
params = []
chi2 = []
for idx,percent in enumerate(percent_bootstrapped):
    if (error[idx]>1e-14).all():
        results = call_model((percent*n).T,time,error[idx].T)
        params.append(results.params)
        chi2.append(results.chisqr)

def plot_values(params,data,time):
    x0 = data[:,0]
    data = data[:,1:]
    M = np.asarray([[params['k1'],0,0],
                    [params['d1'],params['k1']+params['k2'],0],
                    [0,params['d2'],params['k3']]])
    for i in range(np.size(M,axis=0)):
        M[i,i] = 2*M[i,i]-np.sum(M[:,i])
    y = np.zeros((np.size(time),np.size(data,axis=0)))
    for idx,t in enumerate(time):
        y[idx,:] = np.dot(expm(M*t),x0)
    return y

pred_t = np.linspace(0,7,40)
pred_y = plot_values(a0_result.params,(percent_a0.T*n)[[2,3,1],:],pred_t)
act_y = (percent_a0.T*n)[[2,3,1],:]
plt.figure()
plt.plot(time,act_y.T,'o')
plt.plot(pred_t,pred_y,'-')
