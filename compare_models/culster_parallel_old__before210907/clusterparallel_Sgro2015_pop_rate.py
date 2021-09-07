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

from Population_parallel_tools import all_indices, get_n_workers
import matplotlib
# need to have the following to save images on 
#if 'DISPLAY' not in os.environ:
#	matplotlib.use('Agg')
import matplotlib.pyplot as plt
# from Sgro2015_agent_and_pop_FUN import Sgro2015_pop
from Sgro2015_agent_and_pop_FUN import Sgro2015_pop_dir_cpl
# Normalization parameters
from NormParam import *
#%% Sgro 2015
title_font_size = 26
sublabel_font_size = 22
trace_width=3
tick_font_size=15
label_font_size=14

#set parameters
e=0.1; tauA=0.09; tauR=tauA/e; sigma = 0.15 # noise strength
N = 100 # number of cells in the population
c0 = 1.2 #  2.0161 #
g = 0.5 # 1.6129 #
#create logspace
j_arr=np.array([0.5]) # j_arr = np.linspace(0, 1, num=21) # 21 
rho_arr= np.array([1e-3]) #np.logspace(-5, -3, num=3)# rho_arr = np.logspace(-5.5, -3, num=26) # 26 
# Pack parameters together
SgroPopParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':c0,'sigma':sigma,'N':N,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}

alphafval = 2 # cAMP input flow
time_separation = 0 # Whether there is time separation assumption

dt=0.005 ; t_tot=25*Nt_Sgro ; t=list(np.arange(0,t_tot,dt))
nSteps = len(t)
dir_cpl = 1 #coupling between cells and media

def calc_updates_Sgro(j_arr, rho_arr, SgroPopParam, nSteps, index):
    ''' index: (j_index,rho_index) '''
    SgroPopParam['j'] = j_arr[index[0]]
    SgroPopParam['rho'] = rho_arr[index[1]]
    
    A_trace_orig=np.zeros((N,nSteps))
    R_trace_orig=np.zeros((N,nSteps))
    cAMPext_trace = np.zeros((nSteps,1))
    
#    A0 = np.random.uniform(-1.5,2,N) # + np.random.uniform(-sigma,sigma,N); 
#    R0 = np.random.uniform(-0.5,2,N) # + np.random.uniform(-sigma,sigma,N)
    A0 = -1.5*np.ones(N) # np.random.uniform(-1.5,2,N) # 
    R0 = -0.5*np.ones(N) # np.random.uniform(-0.5,2,N) #
    
    cAMPext0 = 0
    A_trace_orig[:,0]=A0 
    R_trace_orig[:,0]=R0
    cAMPext_trace[0] = cAMPext0

    Sgro_pop=Sgro2015_pop_dir_cpl(A0,R0,cAMPext0, SgroPopParam)
    for i in range(nSteps-2):
        A_next,R_next, cAMPext_next = Sgro_pop.update(dt, time_separation, alphafval,dir_cpl)
        A_trace_orig[:,i+1] = A_next
        R_trace_orig[:,i+1] = R_next
        cAMPext_trace[i+1] = cAMPext_next
    # Traces
    A_trace_offset=1.5
    A_trace_plot=(A_trace_orig+A_trace_offset)/SgroPopParam['Na'];
    A_trace_mean = np.mean(A_trace_plot,axis = 0)   
    R_trace_plot = R_trace_orig
    R_trace_mean = np.mean(R_trace_plot,axis = 0)

    later_portion = 0.2 # start count peaks after this X total simulation time
    A_trace_later = A_trace_plot[:,math.floor(nSteps * later_portion):]
    A_trace_mean_later=A_trace_mean[math.floor(nSteps * later_portion):] # the later part of trace
    PkPos, PkProperties = find_peaks(A_trace_mean_later, prominence=(0.4,2))
    cAMPext_trace_later = cAMPext_trace[math.floor(nSteps * later_portion):] 
    t_plot = np.array(t)/Nt_Sgro
    t_plot_later = t_plot[math.floor(nSteps * later_portion):]

    if len(PkPos) == 0:
        firing_rate = 0; height = 0
    else: 
        firing_rate = len(PkPos)/(t_tot/Nt_Sgro*(1-later_portion))
        height = np.mean(PkProperties["prominences"])
    
    # check simulation traces
    fig = plt.figure(figsize=(7,4.5)); grid = plt.GridSpec(1, 1,hspace= 0.3)
    ax1= fig.add_subplot(grid[0, 0])
    # single cell traces
    cell_idx = [0,10,20,30,40,50,60,70,80]
    for i in cell_idx:
#        ax1.plot(t_plot_later,A_trace_later[i,:], color = 'lightgreen', linewidth = 2)
        ax1.plot(t_plot,A_trace_plot[i,:], color = 'lightgreen', linewidth = 2)
        ax1.plot(t_plot,R_trace_plot[i,:], color = 'gold', linewidth = 2)
        
    
#    ax1.plot(t_plot_later,A_trace_mean_later, color ='darkgreen',linewidth = 3,
#             label = 'pop mean cAMPi')
    ax1.plot(t_plot,A_trace_mean, color = 'darkgreen', linewidth = 2)
    ax1.plot(t_plot,R_trace_mean, color = 'orange', linewidth = 2)
    ax2 = ax1.twinx()
    ax2.plot(t_plot,cAMPext_trace, color ='k',linewidth=3,
             label = 'cAMPe')
    ax2.set_ylabel(r'$cAMP_{e}$, A.U.',color = 'k',fontsize=label_font_size) #  fontsize=label_font_size-3
    ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    ax1.set_xlim([0,15])
    ax1.set_xlabel('Time, A.U.')
    ax1.set_ylabel(r'$cAMP_{i}$, A.U.',color = 'g',fontsize=label_font_size)
    ax1.set_title('dilutio rate= '+ '{:#.3n}'.format(np.float64(j_arr[index[0]])) +
    ', density= '+'{:#.3n}'.format(np.float64(rho_arr[index[1]])) + 
        ', FR = '+'{:#.3n}'.format(np.float64(firing_rate)) + 
            ', coupling='+str(dir_cpl)+'not random initcond')
    leg = ax1.legend()
                
    # return index,firing_rate, height
    return index,firing_rate, height, A_trace_plot

# To run the simulation get all indices for gamma and rho
indices = all_indices(j_arr, rho_arr)

# Create a partial function call to calc_updates
par_calc_updates_Sgro = functools.partial(calc_updates_Sgro, j_arr, rho_arr, SgroPopParam, nSteps)

n_workers = get_n_workers()

# Start up the worker pool and return the results as a list of
# [((gamma_index, rho_index), len(PkPos)),...]
tic = perf_counter()
#results = list(map(par_calc_updates_Sgro,indices))
#
## Run serially or in parallel if possible
#if n_workers == 1:
#    results = list(map(par_calc_updates_Sgro,indices))
#else:
#    with mp.Pool(get_n_workers()) as pool:
#        results = pool.map(par_calc_updates_Sgro,indices)

# use just 1 worker
results = list(map(par_calc_updates_Sgro,indices))
        
toc = perf_counter() 
print('time passed with %s workers: %s' % (get_n_workers(),toc-tic))

#%% assign results to pop_rate matrix
pop_rate_Sgro = np.zeros([j_arr.size, rho_arr.size]) # population firing rate
pop_height_Sgro = np.zeros([j_arr.size, rho_arr.size]) 

for i in results:
    idx = i[0]; firing_rate = i[1]; height = i[2]
    pop_rate_Sgro[idx] = firing_rate
    pop_height_Sgro[idx] = height
np.savez('pop_fire_rate_Sgro_OUT_191025.npz', j_arr = j_arr, rho_arr = rho_arr,pop_rate_Sgro = pop_rate_Sgro, pop_height_Sgro=pop_height_Sgro)
    
# Plot heat map
title_font_size = 20 
label_font_size = 20
tick_font_size = 16

fig3 = plt.figure(figsize=(23, 4.5))
grid = plt.GridSpec(1,3, wspace=0.3, hspace=0.05)

ax1= fig3.add_subplot(grid[0,0])
heatmap = ax1.pcolor(j_arr, rho_arr, pop_rate_Sgro.transpose(), cmap='jet') # cmap='jet'
ax1.set_yscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax1);cbar.ax.tick_params(labelsize = tick_font_size) 
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Sgro&Mehta 2015 pop firing rate', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax2= fig3.add_subplot(grid[0,1])
heatmap = ax2.pcolor(j_arr, rho_arr, pop_height_Sgro.transpose(), cmap='jet') # cmap='jet'
ax1.set_yscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}(\gamma)$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Sgro&Mehta 2015 pop firing height, TS='+str(time_separation)+', noise = '+str(sigma),
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# save image and results to the current folder
plot_name = 'pop_fire_rate_Sgro_191025.png'
plt.savefig(plot_name)
plt.show()

