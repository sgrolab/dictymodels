# -*- coding: utf-8 -*-
"""
Test parallel

Important: Run on the SCC in batch or interactive jobs with at least 2 cores (qsub/qrsh flag:  -pe omp N where N is 2,4,8,16,28,32)
For development on a login node, set the NSLOTS variable when running Spyder or Python:
    
    NSLOTS=2 spyder &
    NSLOTS=2 python test_parallel.py
    
If running on Windows or OSX it'll use all available CPUs.
"""

import numpy as np
from time import perf_counter 
import functools
import multiprocessing as mp
from scipy.signal import find_peaks
import math

from Population_parallel_tools import all_indices, get_n_workers
import matplotlib.pyplot as plt

import scipy.io
import pandas as pd

# needs to be determined
Nt_Gregor = 6
Nt_Sgro = 27
Nt_Goldbeter = 7
Nt_Maeda = 3.57
Nt_Kamino = 6
#%%  Kamino 2017
from Kamino2017_agent_and_pop_FUN import Kamino2017_pop
# Set some parameters
tau=1.5  
n=2  
K=4
kt=2
delta=0.01

#gamma_arr=np.logspace(0, 2.0, num=2) # 21
#rho_arr=np.logspace(0, 2.0, num=3) # 26
gamma_arr=np.logspace(0, 2.0, num=5) # 21
rho_arr=np.array([10])

# Pack parameters together
Param={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta}


# Initialize oscillation matrix
pop_rate_Kamino = np.zeros([gamma_arr.size, rho_arr.size]) # population firing rate


dt=0.001  
t_tot=10*Nt_Kamino 
# Number of time steps
t=np.arange(0,t_tot,dt)
nSteps = len(t)
t_plot_Kamino = np.array(t)
z0_influx = 0

def calc_updates_Kamino(gamma_arr, rho_arr, Param, nSteps, index):
    ''' index: (gamma_index,rho_index) '''
    Param['gamma'] = gamma_arr[index[0]]
    Param['rho'] = rho_arr[index[1]]
    
    x0=0.01
    y0=0.08
    z0=0
    x_trace=[x0]
    y_trace=[y0]
    z_trace=[z0]

    Kamino_pop=Kamino2017_pop([x0,y0,z0],Param)
    for i in range(len(t)-1):      
        x_next,y_next,z_next=Kamino_pop.update(z0_influx,dt)
        x_trace.append(x_next)
        y_trace.append(y_next)
        z_trace.append(z_next)
            
    y_trace=np.array(y_trace) # convert list to array
    later_portion = 0.2 # start count peaks after this X total simulation time
    y_trace_later=y_trace[math.floor(nSteps * later_portion):] # the later part of trace
        
    # check simulation traces
    fig = plt.figure()
    plt.plot(t_plot_Kamino,x_trace,t_plot_Kamino,y_trace,t_plot_Kamino,z_trace)
    plt.xlabel('Time')
    plt.ylabel('x,y,z')
    plt.title('Fig5D with gamma= '+str(gamma_arr[index[0]])+' rho= '+str(rho_arr[index[1]]))
    plt.gca().legend(('x','y','z'))
    plt.show()

    PkPos, PkProperties = find_peaks(y_trace_later, prominence=(0.02,100))
#    # Check find_peaks
#    fig = plt.figure()
#    plt.plot(y_trace_later)
#    plt.plot(PkPos, y_trace_later[PkPos], "x")
    
    if len(PkPos) == 0:
        firing_rate = 0
    else: 
        firing_rate = len(PkPos)/(t_tot/Nt_Kamino*(1-later_portion))
    return index,firing_rate

# To run the simulation get all indices for gamma and rho
indices = all_indices(gamma_arr, rho_arr)

# Create a partial function call to calc_updates
par_calc_updates_Kamino = functools.partial(calc_updates_Kamino, gamma_arr, rho_arr, Param, nSteps)

n_workers = get_n_workers()

# Start up the worker pool and return the results as a list of
# [((gamma_index, rho_index), len(PkPos)),...]
tic = perf_counter()

# use 1 worker
results = list(map(par_calc_updates_Kamino,indices))

## Run serially or in parallel if possible
#if n_workers == 1:
#    results = list(map(par_calc_updates,indices))
#else:
#    with mp.Pool(get_n_workers()) as pool:
#        results = pool.map(par_calc_updates,indices)

toc = perf_counter() 
print('time passed with %s workers: %s' % (get_n_workers(),toc-tic))

#%% assign results to pop_rate matrix
for i in results:
    idx = i[0]; firing_rate = i[1]
    pop_rate_Kamino[idx] = firing_rate
    
# Plot heat map
title_font_size = 20 
label_font_size = 20
tick_font_size = 16

fig3 = plt.figure(figsize=(23, 4.5))
grid = plt.GridSpec(1,3, wspace=0.3, hspace=0.05)


ax2= fig3.add_subplot(grid[0,1])
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, OscOrNot.transpose(), cmap='jet') # cmap='jet'
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
heatmap = ax2.pcolor(gamma_arr, rho_arr, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
# ax2.set_xscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}(\gamma)$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Kamino & Sawai 2017', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# results_dir = 'C:/Users/ellin/Documents/GitHub/dictymodels/Python agent models/cluster parallel/'
plot_name = 'pop_fire_rate_Kamino2017'
plt.savefig(plot_name)
# plt.savefig(result_dir + plot_name)
np.savez('pop_fire_rate_Kamino_OUT.npz', gamma_arr = gamma_arr, rho_arr = rho_arr,pop_rate_Kamino = pop_rate_Kamino)


#%% Sgro 2015
from Sgro2015_agent_and_pop_FUN import Sgro2015_pop
#set parameters
e=0.1; tauA=0.09; tauR=tauA/e; g=0.5; sigma = 0.15 # noise strength
N = 100 # number of cells in the population

#create logspace
j_arr=np.linspace(0, 1, num=2) # 21
rho_arr=np.logspace(-5.5, -3, num=2) # 26

SgroPopParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':sigma,'N':N,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
# Initialize firing rate  matrix
pop_rate_Sgro = np.zeros([j_arr.size, rho_arr.size]) # population firing rate

dt=0.005 ; t_tot=200*Nt_Sgro ; t=list(np.arange(0,t_tot,dt))
nSteps = len(t)
t_plot_Sgro = np.array(t); t_plot_Sgro = t_plot_Sgro/SgroPopParam['Nt']
# cAMP input flow
alphafval = 2.5
# Whether there is time separation assumption
time_separation = 0

#update function

def calc_updates_Sgro(j_arr, rho_arr, SgroPopParam, nSteps, index):
    ''' index: (j_index,rho_index) '''
    SgroPopParam['j'] = j_arr[index[0]]
    SgroPopParam['rho'] =  rho_arr[index[1]]
    
    A_trace_orig=np.zeros((N,nSteps))
    R_trace_orig=np.zeros((N,nSteps))
    cAMPext_trace = np.zeros((nSteps,1))
    
    A0=-1.5*np.ones(N) # + np.random.uniform(-sigma,sigma,N); 
    R0=-0.5*np.ones(N) # + np.random.uniform(-sigma,sigma,N); 
    cAMPext0 = 0

    A_trace_orig[:,0]=A0 
    R_trace_orig[:,0]=R0
    cAMPext_trace[0] = cAMPext0
    Sgro_pop=Sgro2015_pop(A0,R0,cAMPext0, SgroPopParam)
    
    for i in range(len(t)-1):
        A_next,R_next, cAMPext_next = Sgro_pop.update(dt, time_separation, alphafval)
        A_trace_orig[:,i] = A_next
        R_trace_orig[:,i] = R_next
        cAMPext_trace[i] = cAMPext_next
    # Traces
    A_trace_offset=1.5
    A_trace_plot=(A_trace_orig+A_trace_offset)/SgroPopParam['Na'];
    A_trace_mean = np.mean(A_trace_plot,axis = 0)

    later_portion = 0.2 # start count peaks after this X total simulation time
    A_trace_mean_later=A_trace_mean[math.floor(nSteps * later_portion):] # the later part of trace
        
#    # check simulation traces
#    fig = plt.figure()
#    plt.plot(t_plot_Sgro,A_trace_mean)
#    plt.xlabel('Time')
#    plt.ylabel('A')
#    plt.title('j= '+str(j_arr[index[0]])+' rho= '+str(rho_arr[index[1]]))
#    plt.show()
#
    PkPos, PkProperties = find_peaks(A_trace_mean_later, prominence=(0.5,2))
#    # Check find_peaks
#    fig1 = plt.figure()
#    plt.plot(A_trace_mean_later)
#    plt.plot(PkPos, A_trace_mean_later[PkPos], "x")
    
    if len(PkPos) == 0:
        firing_rate = 0
    else: 
        firing_rate = len(PkPos)/(t_tot/Nt_Sgro*(1-later_portion))
    return index,firing_rate


# To run the simulation get all indices for gamma and rho
indices = all_indices(j_arr, rho_arr)

# Create a partial function call to calc_updates
par_calc_updates_Sgro = functools.partial(calc_updates_Sgro, j_arr, rho_arr, SgroPopParam, nSteps)

n_workers = get_n_workers()

# Start up the worker pool and return the results as a list of
# [((j_index, rho_index), len(PkPos)),...]
tic = perf_counter()

# use 1 worker
results = list(map(par_calc_updates_Sgro,indices))

## Run serially or in parallel if possible
#if n_workers == 1:
#    results = list(map(par_calc_updates_Sgro,indices))
#else:
#    with mp.Pool(get_n_workers()) as pool:
#        results = pool.map(par_calc_updates_Sgro,indices)

toc = perf_counter() 
print('time passed with %s workers: %s' % (get_n_workers(),toc-tic))

# assign results to pop_rate matrix
for i in results:
    idx = i[0]; firing_rate = i[1]
    pop_rate_Sgro[idx] = firing_rate
    
# Plot heat map
title_font_size = 20 
label_font_size = 20
tick_font_size = 16

fig3 = plt.figure(figsize=(23, 4.5))
grid = plt.GridSpec(1,3, wspace=0.3, hspace=0.05)


ax2= fig3.add_subplot(grid[0,1])
heatmap = ax2.pcolor(j_arr, rho_arr, pop_rate_Sgro.transpose(), cmap='jet') # cmap='jet'
# ax2.set_xscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}j$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Sgro 2015', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# results_dir = 'C:/Users/ellin/Documents/GitHub/dictymodels/Python agent models/cluster parallel/'
plot_name = 'pop_fire_rate_Sgro2015'
plt.savefig(plot_name)
# plt.savefig(result_dir + plot_name)

np.savez('pop_fire_rate_Sgro_OUT.npz', j_arr = j_arr, rho_arr = rho_arr,pop_rate_Sgro = pop_rate_Sgro)

#%% Golbeter 1987
from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_pop_3var
# set parameters
k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10; L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01; theta=0.01
e=  1 # 0.108 # compared to 1
q=4000
sig= 0.6 # 0.57 
v=12; k= 4 # k prime in the paper
ki=1.7 # 0.958 
kt=0.9

#kc=5.4 # 3.58 #j
#h=5 #rho

kc_arr = np.linspace(0.1, 25, num=50)
h_arr = np.linspace(0.1, 25, num=50)

Goldbeter3PopParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':0,'h':0}

#update function
dt=0.001; t_tot=200*Nt_Goldbeter; t=list(np.arange(0,t_tot,dt))
nSteps = len(t)
t_plot_Goldbeter = np.array(t); t_plot_Goldbeter = t_plot_Goldbeter/Nt_Goldbeter
signal_input = 0

def calc_updates_Goldbeter(kc_arr, h_arr, Goldbeter3PopParam, nSteps, index):
    ''' index: (kc_index,h_index) '''
    Goldbeter3PopParam['kc'] = kc_arr[index[0]]
    Goldbeter3PopParam['h'] =  h_arr[index[1]]
    

    p0=0.8; a0=3; b0=0.9; g0=0
    p_trace=[p0]; b_trace=[b0]; g_trace=[g0]

    Goldbeter3_pop=Goldbeter1987_pop_3var(0,[p0,a0,b0,g0],Goldbeter3PopParam)
    for i in range(len(t)-1):
        p_next,b_next,g_next= Goldbeter3_pop.update(dt,a0,signal_input)
        p_trace.append(p_next)
        b_trace.append(b_next)
        g_trace.append(g_next)
        
   
    # Convert into np array
    b_trace = np.array(b_trace)
    later_portion = 0.2 # start count peaks after this X total simulation time
    b_trace_later=b_trace[math.floor(nSteps * later_portion):] # the later part of trace
    # check simulation traces
#    fig = plt.figure()
#    plt.plot(t_plot_Goldbeter,b_trace)
#    plt.xlabel('Time')
#    plt.ylabel('b')
#    plt.title('kc= '+str(kc_arr[index[0]])+' h= '+str(h_arr[index[1]]))
#    plt.show()

    PkPos, PkProperties = find_peaks(b_trace_later, prominence=(5,1000))
#    # Check find_peaks
#    fig = plt.figure()
#    plt.plot(b_trace_later)
#    plt.plot(PkPos, b_trace_later[PkPos], "x")
    
    if len(PkPos) == 0:
        firing_rate = 0
    else: 
        firing_rate = len(PkPos)/(t_tot/Nt_Goldbeter*(1-later_portion))
    return index,firing_rate

# To run the simulation get all indices for gamma and rho
indices = all_indices(kc_arr, h_arr)

# Create a partial function call to calc_updates
par_calc_updates_Goldbeter = functools.partial(calc_updates_Goldbeter, kc_arr, h_arr, Goldbeter3PopParam, nSteps)

n_workers = get_n_workers()

# Start up the worker pool and return the results as a list of
# [((gamma_index, rho_index), len(PkPos)),...]
tic = perf_counter()

# use 1 worker
#results = list(map(par_calc_updates_Goldbeter,indices))

# Run serially or in parallel if possible
if n_workers == 1:
    results = list(map(par_calc_updates_Goldbeter,indices))
else:
    with mp.Pool(get_n_workers()) as pool:
        results = pool.map(par_calc_updates_Goldbeter,indices)

toc = perf_counter() 
print('time passed with %s workers: %s' % (get_n_workers(),toc-tic))


#plot heat map
# assign results to pop_rate matrix

pop_rate_Goldbeter = np.zeros([kc_arr.size, h_arr.size])
for i in results:
    idx = i[0]; firing_rate = i[1]
    pop_rate_Goldbeter[idx] = firing_rate
    
# Plot heat map
title_font_size = 20 
label_font_size = 20
tick_font_size = 16

fig3 = plt.figure(figsize=(6, 4.5))
grid = plt.GridSpec(1,3, wspace=0.3, hspace=0.05)


ax2= fig3.add_subplot(grid[0,1])
heatmap = ax2.pcolor(kc_arr, h_arr, pop_rate_Goldbeter.transpose(), cmap='jet') # cmap='jet'
# ax2.set_xscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}j$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Goldbeter 1987', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
plt.show()

# save image and results to the current folder
# results_dir = 'project/sgrolab/cqhuyan/test_run_time/Gold/'
plot_name = '~/project/sgrolab/cqhuyan/test_run_time/Gold/pop_fire_rate_Goldbeter1987.png'
plt.savefig(plot_name)
# plt.savefig(result_dir + plot_name)
np.savez('pop_fire_rate_Goldbeter_OUT.npz', kc_arr = kc_arr, h_arr = h_arr,pop_rate_Goldbeter = pop_rate_Goldbeter)

#%%  Maeda & Loomis 2004
from MaedaLoomis2004_agent_and_pop_FUN import MaedaLoomis2004_pop
# Set some parameters

# parameters from Maeda & Loomis 2004 paper
k1=2.0; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=1.0; k8=1.3; k9=0.3; k10=0.8
k11=0.7; k12=4.9; k13=23; k14=4.5

# Pack parameters together
MaedaAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}

#gamma_arr=np.logspace(0, 1.6, num=20) # 21
# rho_arr=np.logspace(-0.1, 2.4, num=20) # 26
gamma_arr = np.array([3.3])
rho_arr=np.array([2])

# Initialize oscillation matrix
pop_rate_Maeda = np.zeros([gamma_arr.size, rho_arr.size]) # population firing rate
pop_height_Maeda = np.zeros([gamma_arr.size, rho_arr.size]) # population firing rate

dt=0.001 # 0.00005
t_tot=50*Nt_Maeda
# Number of time steps
t=np.arange(0,t_tot,dt)
nSteps = len(t)
t_plot_Maeda = np.array(t)
signal_input = 0

def calc_updates_Maeda(gamma_arr, rho_arr, nSteps, index):
    ''' index: (gamma_index,rho_index) '''
    gamma = gamma_arr[index[0]]
    rho = rho_arr[index[1]]
    
    ACA0=0.1; PKA0=0.1; ERK20=0.1; RegA0=0.1; cAMPi0=0.01; 
    cAMPe0=0.1; CAR10=0.1
    state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]
    
    ACA_trace=[ACA0]; PKA_trace=[PKA0]; ERK2_trace=[ERK20]
    RegA_trace=[RegA0]; cAMPi_trace=[cAMPi0]; cAMPe_trace=[cAMPe0]
    CAR1_trace=[CAR10]

    MaedaLoomis_pop=MaedaLoomis2004_pop([1,1],state0,MaedaAgentParam)
    for i in range(len(t)-1):      
        ACA_next,PKA_next,ERK2_next,RegA_next,\
        cAMPi_next,cAMPe_next,CAR1_next=MaedaLoomis_pop.update(dt,signal_input,rho,gamma)
        ACA_trace.append(ACA_next)
        PKA_trace.append(PKA_next)
        ERK2_trace.append(ERK2_next)
        RegA_trace.append(RegA_next)
        cAMPi_trace.append(cAMPi_next)
        cAMPe_trace.append(cAMPe_next)
        CAR1_trace.append(CAR1_next)
                    
    cAMPi_trace=np.array(cAMPi_trace) # convert list to array
    later_portion = 0.2 # start count peaks after this X total simulation time
    cAMPi_trace_later=cAMPi_trace[math.floor(nSteps * later_portion):] # the later part of trace
        
#    # check simulation traces
#    fig = plt.figure()
#    plt.plot(t_plot_Kamino,y_trace,t_plot_Kamino,z_trace)
#    plt.xlabel('Time')
#    plt.ylabel('x,y,z')
#    plt.title('Fig5D with gamma= '+str(gamma_arr[index[0]])+' rho= '+str(rho_arr[index[1]]))
#    plt.gca().legend(('y','z'))
#    plt.show()

    PkPos, PkProperties = find_peaks(cAMPi_trace_later, prominence=(0.01,10000))
    
#    # Check find_peaks
#    fig = plt.figure(figsize=(5, 5))
#    plt.plot(cAMPi_trace_later)
#    plt.plot(PkPos, cAMPi_trace_later[PkPos], "x")
#    plt.title('gamma= '+str(gamma_arr[index[0]])+' rho= '+str(rho_arr[index[1]]))
#    
    if len(PkPos) == 0:
        firing_rate = 0
        height = 0
    else: 
        firing_rate = len(PkPos)/(t_tot/Nt_Maeda*(1-later_portion))
        height = np.mean(PkProperties["prominences"])
    return index,firing_rate, height,cAMPi_trace, PkProperties

# To run the simulation get all indices for gamma and rho
indices = all_indices(gamma_arr, rho_arr)

# Create a partial function call to calc_updates
par_calc_updates_Maeda = functools.partial(calc_updates_Maeda, gamma_arr, rho_arr, nSteps)

n_workers = get_n_workers()

# Start up the worker pool and return the results as a list of
# [((gamma_index, rho_index), len(PkPos)),...]
tic = perf_counter()

# use 1 worker
results = list(map(par_calc_updates_Maeda,indices))

## Run serially or in parallel if possible
#if n_workers == 1:
#    results = list(map(par_calc_updates,indices))
#else:
#    with mp.Pool(get_n_workers()) as pool:
#        results = pool.map(par_calc_updates,indices)

toc = perf_counter() 
print('time passed with %s workers: %s' % (get_n_workers(),toc-tic))

#%% assign results to pop_rate matrix
for i in results:
    idx = i[0]; firing_rate = i[1];height = i[2]
    pop_rate_Maeda[idx] = firing_rate
    pop_height_Maeda[idx] = height 

# Plot heat map
title_font_size = 20 
label_font_size = 20
tick_font_size = 16

fig3 = plt.figure(figsize=(23, 4.5))
grid = plt.GridSpec(1,3, wspace=0.3, hspace=0.05)


ax1= fig3.add_subplot(grid[0,0])
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, OscOrNot.transpose(), cmap='jet') # cmap='jet'
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
heatmap = ax1.pcolor(gamma_arr, rho_arr, pop_rate_Maeda.transpose(), cmap='jet') # cmap='jet'
# ax2.set_xscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax1);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}(\gamma)$',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Maeda&Loomis2004 pop firing rate', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax2= fig3.add_subplot(grid[0,1])
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, OscOrNot.transpose(), cmap='jet') # cmap='jet'
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
heatmap = ax2.pcolor(gamma_arr, rho_arr, pop_height_Maeda.transpose(), cmap='jet') # cmap='jet'
# ax2.set_xscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}(\gamma)$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Maeda&Loomis2004 pop firing height', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})


# results_dir = 'C:/Users/ellin/Documents/GitHub/dictymodels/Python agent models/cluster parallel/'
plot_name = 'pop_fire_rate_Maeda2017'
plt.savefig(plot_name)
# plt.savefig(result_dir + plot_name)
np.savez('pop_fire_rate_Kamino_OUT.npz', gamma_arr = gamma_arr, rho_arr = rho_arr,pop_rate_Kamino = pop_rate_Kamino)

#%% Simulation outputs 
# Sgro 2015 from Matlab outputs--> xlsx (simulation results from matlab)
OUT_path = 'C:/Users/ellin/Documents/GitHub/dictymodels/Python agent models/pop_fire_rate_OUT/'
pop_rate_path = OUT_path + r'Sgro_pop_rate.xlsx'
pop_rate = pd.read_excel(pop_rate_path)
pop_rate_Sgro=pop_rate.as_matrix()

rho_arr_Sgro = np.logspace(-5.5,-3,num=26)
j_arr_Sgro = np.linspace(0,1,num=21)

Sgro_low_noise_OUT = np.load(OUT_path +'pop_fire_rate_Sgro_OUT_191027_same_init_cond_ttot25_sigma0.1.npz')
j_arr_Sgro_low_noise =  Sgro_low_noise_OUT['j_arr']
rho_arr_Sgro_low_noise =  Sgro_low_noise_OUT['rho_arr']
pop_rate_Sgro_low_noise =  Sgro_low_noise_OUT['pop_rate_Sgro']


# Gregor 2010 from Matlab outputs--> xlsx
pop_rate_path = OUT_path + r'Gregor_pop_rate.xlsx'
pop_rate_Gregor = pd.read_excel(pop_rate_path)
pop_rate_mat_Gregor=pop_rate_Gregor.as_matrix()

rho_arr_Gregor = np.logspace(-3.5,1,num=26)
k_arr_Gregor = np.linspace(1,100,num=21)

#Golbeter 1987 from parallel computing outputs
#Goldbeter_OUT = np.load(OUT_path +'pop_fire_rate_Goldbeter_OUT_191004.npz')
Goldbeter_OUT = np.load(OUT_path +'pop_fire_rate_Goldbeter_OUT_191122_3.npz')
kc_arr_Goldbeter =  Goldbeter_OUT['kc_arr']
h_arr_Goldbeter =  Goldbeter_OUT['one_over_h_arr']
pop_rate_Goldbeter = Goldbeter_OUT['pop_rate_Goldbeter']
pop_height_Goldbeter = Goldbeter_OUT['pop_height_Goldbeter']

#Maeda Loomis 2004 from parallel computing outputs
Maeda_OUT = np.load(OUT_path +'pop_fire_rate_Maeda_OUT_191004.npz')
rho_arr_Maeda =  Maeda_OUT['rho_arr']
gamma_arr_Maeda =  Maeda_OUT['gamma_arr']
pop_rate_Maeda = Maeda_OUT['pop_rate_Maeda']

# Kamino 2017 from parallel computing outputs
Kamino_OUT = np.load(OUT_path +'pop_fire_rate_Kamino_OUT_191005.npz')
rho_arr_Kamino =  Kamino_OUT['rho_arr']
gamma_arr_Kamino =  Kamino_OUT['gamma_arr']
pop_rate_Kamino = Kamino_OUT['pop_rate_Kamino']
#%% Experimental data, needs to be filled!!!
PopRateExp = np.array([[0.129660088,0.171190476,0.152777778,0.171016069,0.132986111,0.149074074,0.144444444,0.177170868,0.162719633,0.178571429],
[0.065434419,0.131937322,0.131658497,0.156811146,0.147619048,0.162799043,0.147368421,0.171820616,0.11875,0.1],
[0.075416667,0.152597403,0.138703704,0.161497326,0.183333333,0.157142857,0.132373581,0.139815266,0.117989418,0.103125],
[0.057236842,0.164138177,0.13287037,0.165909091,0.126666667,0.166827485,0.086222222,0.065789474,0.047584541,0.055555556],
[0.019047619,0.103693161,0.122423896,0.100198413,0.106683375,0.091592443,0.034122807,0.042105263,0.031818182,0],
[0.014285714,0.026428571,0.031578947,0.049938272,0.038095238,0.00625,0.010526316,0,0,0],
[0,0.02375,0.022222222,0.050595238,0,0,0,0,0,0]])

JExp = np.linspace(1,10, num=10) # np.array([0,1,2,4,6,8,10,15,16,20])
RhoExp = np.linspace(1,7, num=7) #  np.array([0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125])

# rho_exp = Sgro2015Fig5Data[0,1:10]
#%% Plot all models
#abcd_font_size = 28
#title_font_size = 25
#label_font_size = 25
#tick_font_size = 20
title_font_size = 18
label_font_size = 18
sublabel_font_size = 14
tick_font_size = 14
trace_width = 2
abcd_font_size = 28

mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']

fig3 = plt.figure(figsize=(24, 8))
grid = plt.GridSpec(2,4, wspace=0.8, hspace=0.5)

ax0= fig3.add_subplot(grid[0,0])
ax0.set_xticks([0,1,2,3,4,5,6,7,8,9]); 
ax0.set_xticklabels([1,2,3,4,5,6,7,8,9,10],fontsize=tick_font_size)
ax0.set_yticks([0,1,2,3,4,5,6,7]); 
ax0.set_yticklabels(['1/2','1/4','1/8','1/16','1/32','1/64','1/128'],fontsize=tick_font_size)
ax0.set_title('Experiment', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax0.set_xlabel('Flow Rate (mL/min)', size=sublabel_font_size)
ax0.set_ylabel('Cell Density(mML)', size=sublabel_font_size)
heatmap = ax0.imshow(PopRateExp, cmap='jet') # cmap='jet'
x=[3.5,4.5,5.5,7.5,9.5]
[ax0.axvline(_x, color='white',linewidth=trace_width) for _x in x]
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax0,ticks=[0,0.05,0.1,0.15]);
cbar.ax.tick_params(labelsize = tick_font_size) 
cbar.set_label( 'cAMP pulses/min',size=tick_font_size)
#ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.text(-0.2 , 1.1, 'A',
         horizontalalignment='center',verticalalignment='center',
         transform = ax0.transAxes, color = 'b', fontsize=abcd_font_size)

ax1= fig3.add_subplot(grid[0,1],xticks=[0,25,50,75,100])
heatmap = ax1.pcolor(kc_arr_Goldbeter, h_arr_Goldbeter, pop_rate_Goldbeter.transpose(), cmap='jet') # cmap='jet'
# ax1.set_xscale('log');
ax1.set_yscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax1, ticks=[0,0.5,1,1.5])
cbar.ax.tick_params(labelsize = tick_font_size) 
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Martiel 1987', color=mycolors[0],fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax1.text(-0.2 , 1.1, 'B',
         horizontalalignment='center',verticalalignment='center',
         transform = ax1.transAxes, color = 'g', fontsize=abcd_font_size)

ax2= fig3.add_subplot(grid[0,2],xticks=[0,25,50,75,100])
heatmap = ax2.pcolor(gamma_arr_Maeda, rho_arr_Maeda, pop_rate_Maeda.transpose(), cmap='jet') # cmap='jet'
ax2.set_yscale('log')
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Maeda 2004', color=mycolors[1],fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax2.text(-0.2 , 1.1, 'C',
         horizontalalignment='center',verticalalignment='center',
         transform = ax2.transAxes, color = 'g', fontsize=abcd_font_size)

ax3= fig3.add_subplot(grid[0,3], xticks=[0,25,50,75,100])
heatmap = ax3.pcolor(k_arr_Gregor, rho_arr_Gregor, pop_rate_Gregor, cmap='jet') # cmap='jet'
ax3.set_yscale('log')
cbar=fig3.colorbar(heatmap, ax=ax3);cbar.ax.tick_params(labelsize = tick_font_size) 
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Gregor 2010',color=mycolors[2], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax3.text(-0.2 , 1.1, 'D',
         horizontalalignment='center',verticalalignment='center',
         transform = ax3.transAxes, color = 'g', fontsize=abcd_font_size)

# Sgro regular noise (sig = 0.15)
ax4= fig3.add_subplot(grid[1,1],xticks=[0,0.25,0.5,0.75,1])
heatmap = ax4.pcolor(j_arr_Sgro, rho_arr_Sgro, pop_rate_Sgro, cmap='jet') # cmap='jet'
ax4.set_yscale('log'); ax4.set_ylim([10**(-5),10**(-3)]); 
cbar=fig3.colorbar(heatmap, ax=ax4);cbar.ax.tick_params(labelsize = tick_font_size) 
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('Sgro 2015', color=mycolors[5],fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax4.text(-0.2 , 1.1, 'E',
         horizontalalignment='center',verticalalignment='center',
         transform = ax4.transAxes, color = 'g', fontsize=abcd_font_size)

# Sgro low noise  (sig = 0.0)
ax5= fig3.add_subplot(grid[1,2],xticks=[0,0.25,0.5,0.75,1])
heatmap = ax5.pcolor(j_arr_Sgro_low_noise, rho_arr_Sgro_low_noise, pop_rate_Sgro_low_noise.transpose(), cmap='jet') # cmap='jet'
ax5.set_yscale('log'); ax5.set_ylim([10**(-5),10**(-3)]); 
cbar=fig3.colorbar(heatmap, ax=ax5);cbar.ax.tick_params(labelsize = tick_font_size) 
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5.set_title('Sgro 2015 \n w/o noise', color=mycolors[5],fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax5.text(-0.2 , 1.1, 'F',
         horizontalalignment='center',verticalalignment='center',
         transform = ax5.transAxes, color = 'g', fontsize=abcd_font_size)

ax6= fig3.add_subplot(grid[1,3], xticks=[0,25,50,75,100])
heatmap = ax6.pcolor(gamma_arr_Kamino, rho_arr_Kamino, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
ax6.set_yscale('log')
cbar=fig3.colorbar(heatmap, ax=ax6);cbar.ax.tick_params(labelsize = tick_font_size) 
ax6.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax6.set_title('Kamino 2017', color=mycolors[7], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax6.text(-0.2 , 1.1, 'G',
         horizontalalignment='center',verticalalignment='center',
         transform = ax6.transAxes, color = 'g', fontsize=abcd_font_size)

# insets
#ax7= fig3.add_axes([0.77,0.13,0.08,0.16])
ax7= fig3.add_axes([0.805,0.13,0.062,0.152])
heatmap = ax7.pcolor(gamma_arr_Kamino, rho_arr_Kamino,pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
ax7.set_yscale('log');ax7.set_xscale('log');
ax7.set_xticks([]) ; ax7.set_yticks([]) 
ax7.spines['bottom'].set_color('white');ax7.spines['top'].set_color('white')
ax7.spines['left'].set_color('white');ax7.spines['right'].set_color('white')

fig3.text(0.62, 0.02, 'Dilution Rate, A.U.',fontsize=label_font_size, ha='center')
fig3.text(0.29, 0.5, 'Population Density, A.U.',fontsize=label_font_size, va='center', rotation='vertical')

plt.show()


