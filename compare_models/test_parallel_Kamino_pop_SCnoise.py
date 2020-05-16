# -*- coding: utf-8 -*-
"""
Test parallel

Important: Run on the SCC in batch or interactive jobs with at least 2 cores (qsub/qrsh flag:  -pe omp N where N is 2,4,8,16,28,32)
For development on a login node, set the NSLOTS variable when running Spyder or Python:
    
    NSLOTS=2 spyder &
    NSLOTS=2 python test_parallel.py
    
If running on Windows or OSX it'll use all available CPUs.
"""
import os
import numpy as np
from time import perf_counter 
import functools
import multiprocessing as mp
from scipy.signal import find_peaks
import math

from parallel_tools import all_indices, get_n_workers
from Kamino2017_agent_and_pop_FUN import Kamino2017_pop_noiseSC
import matplotlib
if 'DISPLAY' not in os.environ:
	matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Normalization parameters
from NormParam import *
Nh_Kamino_offset = 0.0585
Nh_Sgro_offset = -1.5

#%% Kamino 2017 
# Set some parameters
tau=1.5  
n=2  
K=4
kt=2
delta= 0.01
z0_influx = 0
# noise strength
sigma = 0
# number of cells in a population
N=100

ParamArrLen = 25
gamma_arr=np.logspace(0, 2.0, num = ParamArrLen) # 21
rho_arr=np.logspace(0, 2.0, num = ParamArrLen) # 26

# Pack parameters together
Param={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,'sigma':sigma, 'N':N}

dt=0.001 ; t_tot=50*Nt_Kamino; t=list(np.arange(0,t_tot,dt))
nSteps = len(t)


def calc_updates(gamma_arr, rho_arr, Param, nSteps, index):
    ''' index: (gamma_index,rho_index) '''
    Param['gamma'] = gamma_arr[index[0]]
    Param['rho'] = rho_arr[index[1]]
    x0=0.01*np.ones(N)
    y0=0.08*np.ones(N)
    z0=0
    x_trace=np.zeros((N,nSteps)); x_trace[:,0] = x0
    y_trace=np.zeros((N,nSteps)); y_trace[:,0] = y0
    z_trace=np.zeros((nSteps,1)); z_trace[0] = z0
    Kamino_pop = Kamino2017_pop_noiseSC(x0,y0,z0,Param)
    
    for i in range(nSteps-1):
        x_next,y_next,z_next=Kamino_pop.update(z0_influx,dt)
        x_trace[:,i+1] = x_next
        y_trace[:,i+1] = y_next
        z_trace[i+1] = z_next
            
    later_portion = 0.2 # start count peaks after this X total simulation time
    y_trace_mean = np.mean(y_trace,axis=0)
    y_trace_mean_norm = (y_trace_mean-Nh_Kamino_offset)/Nh_Kamino # height normalization
    y_trace_mean_norm_later=y_trace_mean_norm[math.floor(nSteps * later_portion):] # the later part of trace
    
    # PkPos, PkProperties = find_peaks(y_trace_mean_norm_later, prominence=(0.02,100))
    pop_max = np.amax(y_trace_mean_norm_later);
    pop_min = np.amin(y_trace_mean_norm_later)
    pk_find_thresh = 0.03
    PkPos, PkProperties = find_peaks(y_trace_mean_norm_later, prominence = pk_find_thresh)
    if len(PkPos) == 0:
        firing_rate = 0; height = 0
    else: 
        firing_rate = len(PkPos)/(t_tot/Nt_Kamino*(1-later_portion))
        height = np.mean(PkProperties["prominences"])
    return index,firing_rate, height

# To run the simulation get all indices for gamma and rho
indices = all_indices(gamma_arr, rho_arr)

# Create a partial function call to calc_updates
par_calc_updates = functools.partial(calc_updates, gamma_arr, rho_arr, Param, nSteps)

n_workers = get_n_workers()

# Start up the worker pool and return the results as a list of
# [((gamma_index, rho_index), len(PkPos)),...]
tic = perf_counter()
results = list(map(par_calc_updates,indices))

## Run serially or in parallel if possible
if n_workers == 1:
    results = list(map(par_calc_updates,indices))
else:
    with mp.Pool(get_n_workers()) as pool:
        results = pool.map(par_calc_updates,indices)
toc = perf_counter() 
print('time passed with %s workers: %s' % (get_n_workers(),toc-tic))

pk_find_thresh = 0.03

sname = 'pop_FR_Kamino_200417_delta0.01_dt'+str(dt)+'_noise_'+str(sigma)+'_PrmLen_'+str(ParamArrLen)+'PkFindThr'+str(pk_find_thresh)
#%% assign results to pop_rate matrix
pop_rate_Kamino = np.zeros([gamma_arr.size, rho_arr.size]) # population firing rate
pop_height_Kamino = np.zeros([gamma_arr.size, rho_arr.size]) 

for i in results:
    idx = i[0]; firing_rate = i[1]; height = i[2]
    pop_rate_Kamino[idx] = firing_rate
    pop_height_Kamino[idx] = height
np.savez(sname+'.npz', gamma_arr = gamma_arr, rho_arr = rho_arr,pop_rate_Kamino = pop_rate_Kamino, pop_height_Kamino=pop_height_Kamino)
    
# Plot heat map
title_font_size = 20 
label_font_size = 20
tick_font_size = 16

fig3 = plt.figure(figsize=(23, 4.5))
grid = plt.GridSpec(1,3, wspace=0.3, hspace=0.05)

ax1= fig3.add_subplot(grid[0,0])
heatmap = ax1.pcolor(gamma_arr, rho_arr, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
ax1.set_xscale('log'); 
ax1.set_yscale('log');
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax1);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}j$',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Kamino 2017 pop firing rate\n'+', dt '+str(dt)+',noise '+str(sigma), fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax2= fig3.add_subplot(grid[0,1])
heatmap = ax2.pcolor(gamma_arr, rho_arr, pop_height_Kamino.transpose(), cmap='jet') # cmap='jet'
ax2.set_xscale('log'); ax2.set_yscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}(\gamma)$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Kamino2017 \n pop firing height', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# save image and results to the current folder
plot_name = sname+'.png'
plt.savefig(plot_name)

plt.show()

