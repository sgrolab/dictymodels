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

from parallel_tools import all_indices, get_n_workers
import matplotlib.pyplot as plt

#%%  Kamino 2017
from Kamino2017_agent_and_pop_FUN import Kamino2017_pop
# Set some parameters
tau=1.5  
n=2  
K=4
kt=2
delta=0.01

gamma_space=np.logspace(0, 2.0, num=2) # 21
loggamma_space_Kamino=np.log10(gamma_space)
rho_space=np.logspace(0, 2.0, num=3) # 26
logrho_space_Kamino = np.log10(rho_space)

# Pack parameters together
Param={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta}


# Initialize oscillation phase matrix, based on z trace
PkWdthMean = np.zeros([gamma_space.size, rho_space.size]) # PkWdthMean- mean oscillation time
PkPrmMean = np.zeros([gamma_space.size, rho_space.size]) # PkPrmMean - mean oscillation peak prominence

OscOrNot = np.zeros([gamma_space.size, rho_space.size])# OscOrNot - 1 or 0, oscillatory/ nonoscillatory
pop_rate_Kamino = np.zeros([gamma_space.size, rho_space.size]) # population firing rate


dt=0.001 
t_tot=150
# Number of time steps
t=np.arange(0,t_tot,dt)
t_plot_Kamino = np.array(t)
z0_influx = 0


def calc_updates_Kamino(gamma_arr, rho_arr, Param, nSteps, index):
    ''' index: (gamma_index,rho_index) '''
    Param['gamma'] = gamma_arr[index[0]]
    Param['rho'] = rho_arr[index[1]]

    z0_influx = Param['z0_influx'] 
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
    plt.plot(t_plot_Kamino,y_trace,t_plot_Kamino,z_trace)
    plt.xlabel('Time')
    plt.ylabel('x,y,z')
    plt.title('Fig5D with gamma= '+str(gamma_arr[index[0]])+' rho= '+str(rho_arr[index[1]]))
    plt.gca().legend(('y','z'))
    plt.show()

    PkPos, PkProperties = find_peaks(y_trace_later, prominence=(0.02,100))
    # Check find_peaks
    fig = plt.figure()
    plt.plot(y_trace_later)
    plt.plot(PkPos, y_trace_later[PkPos], "x")
    
    if len(PkPos) == 0:
        firing_rate = 0
    else: 
        firing_rate = len(PkPos)/t_tot
    return index,firing_rate

# To run the simulation get all indices for gamma and rho
indices = all_indices(gamma_space, rho_space)

# Create a partial function call to calc_updates
par_calc_updates_Kamino = functools.partial(calc_updates_Kamino, gamma_space, rho_space, Param, nSteps)

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

# assign results to pop_rate matrix
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
heatmap = ax2.pcolor(gamma_space, logrho_space_Kamino, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
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

#%% Sgro 2015


#%% Plot all models
title_font_size = 20 
label_font_size = 20
tick_font_size = 16

fig3 = plt.figure(figsize=(23, 4.5))
grid = plt.GridSpec(1,3, wspace=0.3, hspace=0.05)

#ax1= fig3.add_subplot(grid[0,0])
#heatmap = ax1.pcolor(k_Gregor, logrho_Gregor, pop_rate_mat_Gregor, cmap='jet') # cmap='jet'
#heatmap.set_clim(0,0.16)
#cbar=fig3.colorbar(heatmap, ax=ax1);cbar.ax.tick_params(labelsize = tick_font_size) 
## ax1.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
## ax1.set_xlabel('Dilution rate J',fontsize=label_font_size)
#ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
#ax1.set_title('Gregor 2010', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax2= fig3.add_subplot(grid[0,1])
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, OscOrNot.transpose(), cmap='jet') # cmap='jet'
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
heatmap = ax2.pcolor(gamma_space, logrho_space_Kamino, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
# ax2.set_xscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}(\gamma)$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Kamino & Sawai 2017', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax3= fig3.add_subplot(grid[0,2])
heatmap = ax3.pcolor(100*j_Sgro, logrho_Sgro, pop_rate_mat_Sgro, cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.6)
cbar=fig3.colorbar(heatmap, ax=ax3);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax3.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax3.set_xlabel('Dilution rate',fontsize=label_font_size)
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Sgro 2015', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# fig3.text(0.5, 0.04,'Dilution rate', ha='center', va='center',fontsize=8)
fig3.text(0.06, 0.5, r'Cell density, $log_{10}(\rho)$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)


plt.show()
