# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 18:26:36 2019

@author: Chuqiao
"""
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
#import matplotlib
#if 'DISPLAY' not in os.environ:
#	matplotlib.use('agg')
import matplotlib.pyplot as plt

# needs to be determined
Nt_Gregor = 6
Nt_Sgro = 27
Nt_Goldbeter = 7
Nt_Maeda = 3.57
Nt_Kamino = 6
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

#gamma_arr=np.logspace(0, 2.0, num=40) # 21
#rho_arr=np.logspace(-0.1, 2.7, num=40) # 26
gamma_arr=np.array([15]) # 21
rho_arr=np.array([5]) # np.linspace(3,5，, 4, num=2) # 26


# Initialize oscillation matrix
pop_rate_Maeda = np.zeros([gamma_arr.size, rho_arr.size]) # population firing rate
pop_height_Maeda = np.zeros([gamma_arr.size, rho_arr.size]) # population firing rate

dt=0.0001 
t_tot=30*Nt_Maeda
# Number of time steps
t=np.arange(0,t_tot,dt)
nSteps = len(t)
t_plot_Maeda = np.array(t)
z0_influx = 0

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
    
    signal_input = 0
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
                    
    cAMPi_trace = np.array(cAMPi_trace) # convert list to array
    ERK2_trace = np.array(ERK2_trace)
    cAMPe_trace = np.array(cAMPe_trace)
    
    later_portion = 0.2 # start count peaks after this X total simulation time
    cAMPi_trace_later=cAMPi_trace[math.floor(nSteps * later_portion):] # the later part of trace   
    cAMPi_trace_later_norm = cAMPi_trace_later/np.amax(cAMPi_trace_later)
    ERK2_trace_later = ERK2_trace[math.floor(nSteps * later_portion):]
    ERK2_trace_later_norm = ERK2_trace_later/np.amax(ERK2_trace_later)
    cAMPe_trace_later=cAMPe_trace[math.floor(nSteps * later_portion):] # the later part of trace   
    cAMPe_trace_later_norm = cAMPe_trace_later/np.amax(cAMPe_trace_later)
    t_plot_Maeda_later = t_plot_Maeda[math.floor(nSteps * later_portion):]
    
    
    PkPos, PkProperties = find_peaks(cAMPi_trace_later, prominence=(0.01,2000))
    
#    # Check find_peaks
#    fig = plt.figure()
#    plt.plot(cAMPi_trace_later)
#    plt.plot(PkPos, cAMPi_trace_later[PkPos], "x")
#    plt.title('gamma= '+str(gamma_arr[index[0]])+' rho= '+str(rho_arr[index[1]]))
    
    if len(PkPos) == 0:
        firing_rate = 0; height = 0
    else: 
        firing_rate = len(PkPos)/(t_tot/Nt_Maeda*(1-later_portion))
        height = np.mean(PkProperties["prominences"])
        
    # check simulation traces
    fig = plt.figure(figsize=(5,2.5)); grid = plt.GridSpec(1, 1,hspace= 0.3)
    ax1= fig.add_subplot(grid[0, 0])
    ax1.plot(t_plot_Maeda_later,cAMPi_trace_later, color = 'g', linewidth = 3, label = 'cAMPi')
    ax1.plot(t_plot_Maeda_later,ERK2_trace_later, color ='b',label = 'ERK2')
    ax1.plot(t_plot_Maeda_later,cAMPe_trace_later, color ='k', label = 'cAMPe')
    # ax1.set_xlim([10,20])
    ax1.set_xlabel('Time')
    ax1.set_title('dilution rate= '+ '{:#.3n}'.format(np.float64(gamma_arr[index[0]])) +
    ', density= '+'{:#.3n}'.format(np.float64(rho_arr[index[1]])) + 
        ', FR = '+'{:#.3n}'.format(np.float64(firing_rate)))
    leg = ax1.legend()
                
#    ax2= fig.add_subplot(grid[0, 1])
#    ax2.plot(t_plot_Maeda_later,cAMPi_trace_later_norm, color = 'g', linewidth = 3)
#    ax2.plot(t_plot_Maeda_later,ERK2_trace_later_norm, color ='b')
#    ax2.plot(t_plot_Maeda_later,cAMPe_trace_later_norm, color ='k')
#    ax2.set_xlim([10,20])
#    ax2.set_xlabel('Time')
#    ax2.set_ylabel('Normalized trace, A.U.')
#    ax2.set_title('Normalized traces, zoomed in')

    plt.show()
    
    return index,firing_rate,height

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
#    results = list(map(par_calc_updates_Maeda,indices))
#else:
#    with mp.Pool(get_n_workers()) as pool:
#        results = pool.map(par_calc_updates_Maeda,indices)

toc = perf_counter() 
print('time passed with %s workers: %s' % (get_n_workers(),toc-tic))

#%% assign results to pop_rate matrix
for i in results:
	idx = i[0]; firing_rate = i[1];height = i[2]
	pop_rate_Maeda[idx] = firing_rate
	pop_height_Maeda[idx] = height 
	
np.savez('pop_fire_rate_Maeda_OUT_191004.npz', gamma_arr = gamma_arr, rho_arr = rho_arr,pop_rate_Maeda = pop_rate_Maeda,pop_height_Maeda = pop_height_Maeda)

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
# plt.show()  Need to comment out plt.show to successfully save image

# save image and results to the current folder
plot_name = 'pop_fire_rate_Maeda2004_191004'
plt.savefig(plot_name)
# plt.savefig(result_dir + plot_name)

