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
import matplotlib
if 'DISPLAY' not in os.environ:
	matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Normalization parameters
from NormParam import *

#%%  Maeda & Loomis 2004
from MaedaLoomis2004_agent_and_pop_FUN import MaedaLoomis2004_pop_SCnoise
# Set some parameters

# parameters from Maeda & Loomis 2004 paper
k1=2.0; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=1.0; k8=1.3; k9=0.3; k10=0.8
k11=0.7; k12=4.9; k13=23; k14=4.5

sigma = 0# noise strength
N = 50  # number of cells

# Pack parameters together
MaedaAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14, 'sigma':sigma, 'N':N}

ParamArrLen = 25 # 40
gamma_arr=np.logspace(0, 2.0, num=ParamArrLen) # 40
rho_arr=np.logspace(-0.1, 2.7, num=ParamArrLen) # 40

# Initialize oscillation matrix
pop_rate_Maeda = np.zeros([gamma_arr.size, rho_arr.size]) # population firing rate
pop_height_Maeda = np.zeros([gamma_arr.size, rho_arr.size]) # population firing rate

dt=0.0001 # 0.00005 
t_tot=60*Nt_Maeda
# Number of time steps
t=np.arange(0,t_tot,dt)
nSteps = len(t)
t_plot_Maeda = np.array(t)/Nt_Maeda
z0_influx = 0

def calc_updates_Maeda(gamma_arr, rho_arr, nSteps, index):
    ''' index: (gamma_index,rho_index) '''
    gamma = gamma_arr[index[0]]
    rho = rho_arr[index[1]]
    
    # ACA0=0.1*np.ones((N,1)); PKA0=0.1*np.ones((N,1)); ERK20=0.1*np.ones((N,1)); 
    # RegA0=0.1*np.ones((N,1)); cAMPi0=0.01*np.ones((N,1)); cAMPe0=0.1; CAR10=0.1*np.ones((N,1))
    ACA0=0.1*np.ones(N); PKA0=0.1*np.ones(N); ERK20=0.1*np.ones(N); 
    RegA0=0.1*np.ones(N); cAMPi0=0.01*np.ones(N); cAMPe0=0.1; CAR10=0.1*np.ones(N)
    
    ACA_trace=np.zeros((N,nSteps)); ACA_trace[:,0] = ACA0
    PKA_trace=np.zeros((N,nSteps)); PKA_trace[:,0] = PKA0
    ERK2_trace=np.zeros((N,nSteps)); ERK2_trace[:,0] = ERK20
    RegA_trace= np.zeros((N,nSteps)); RegA_trace[:,0] = RegA0
    cAMPi_trace= np.zeros((N,nSteps)); cAMPi_trace[:,0] = cAMPi0
    cAMPe_trace=np.zeros((nSteps,1));cAMPe_trace[0] = cAMPe0;
    CAR1_trace= np.zeros((N,nSteps)); CAR1_trace[:,0] = CAR10
    
    signal_input = 0
    MaedaLoomis_pop=MaedaLoomis2004_pop_SCnoise([1,1],ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10, MaedaAgentParam)
    for i in range(len(t)-1):      
        ACA_next,PKA_next,ERK2_next,RegA_next,\
        cAMPi_next,cAMPe_next,CAR1_next=MaedaLoomis_pop.update(dt,signal_input,rho,gamma)
        ACA_trace[:,i+1]= ACA_next
        PKA_trace[:,i+1]= PKA_next
        ERK2_trace[:,i+1]= ERK2_next
        RegA_trace[:,i+1]= RegA_next
        cAMPi_trace[:,i+1]= cAMPi_next
        cAMPe_trace[i+1]= cAMPe_next
        CAR1_trace[:,i+1]= CAR1_next
                    
    cAMPi_trace_norm = cAMPi_trace / Nh_Maeda
    cAMPi_trace_norm_mean = np.mean(cAMPi_trace_norm,axis = 0) # population mean
    later_portion = 0.5 # start count peaks after this X total simulation time
    cAMPi_trace_norm_mean_later=cAMPi_trace_norm_mean[math.floor(nSteps * later_portion):] # the later part of trace
    
    pop_max = np.amax(cAMPi_trace_norm_mean_later)
    pop_min = np.amin(cAMPi_trace_norm_mean_later)
    pk_find_prm = 0.7
    PkPos, PkProperties = find_peaks(cAMPi_trace_norm_mean_later,prominence=((pop_max-pop_min)*pk_find_prm,pop_max))
#    # check simulation traces
#    fig = plt.figure()
#    plt.plot(t_plot_Kamino,y_trace,t_plot_Kamino,z_trace)
#    plt.xlabel('Time')
#    plt.ylabel('x,y,z')
#    plt.title('Fig5D with gamma= '+str(gamma_arr[index[0]])+' rho= '+str(rho_arr[index[1]]))
#    plt.gca().legend(('y','z'))
#    plt.show()
	
    # PkPos, PkProperties = find_peaks(cAMPi_trace_mean_norm_later, prominence=(0.01/Nh_Maeda,2000/Nh_Maeda))
    
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
    return index,firing_rate,height

# To run the simulation get all indices for gamma and rho
indices = all_indices(gamma_arr, rho_arr)

# Create a partial function call to calc_updates
par_calc_updates_Maeda = functools.partial(calc_updates_Maeda, gamma_arr, rho_arr, nSteps)

n_workers = get_n_workers()

# Start up the worker pool and return the results as a list of
# [((gamma_index, rho_index), len(PkPos)),...]
tic = perf_counter()

## use 1 worker
#results = list(map(par_calc_updates_Maeda,indices))

# Run serially or in parallel if possible
if n_workers == 1:
    results = list(map(par_calc_updates_Maeda,indices))
else:
    with mp.Pool(get_n_workers()) as pool:
        results = pool.map(par_calc_updates_Maeda,indices)

toc = perf_counter() 
print('time passed with %s workers: %s' % (get_n_workers(),toc-tic))

#%% assign results to pop_rate matrix

save_name = 'pop_FR_Maeda_200331_hnorm_dt'+str(dt)+'_noise'+str(sigma)+'ParamLen'+str(ParamArrLen)

for i in results:
	idx = i[0]; firing_rate = i[1];height = i[2]
	pop_rate_Maeda[idx] = firing_rate
	pop_height_Maeda[idx] = height 
	
np.savez(save_name+'p.npz', gamma_arr = gamma_arr, rho_arr = rho_arr,pop_rate_Maeda = pop_rate_Maeda,pop_height_Maeda = pop_height_Maeda)

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
ax1.set_yscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax1);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}(\gamma)$',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Maeda 2004 pop firing rate'+', dt '+str(dt)+',noise '+str(sigma), fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax2= fig3.add_subplot(grid[0,1])
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, OscOrNot.transpose(), cmap='jet') # cmap='jet'
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
heatmap = ax2.pcolor(gamma_arr, rho_arr, pop_height_Maeda.transpose(), cmap='jet') # cmap='jet'
ax2.set_yscale('log')
# ax2.set_xscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}(\gamma)$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Maeda2004 \n pop firing height', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})


# save image and results to the current folder
plot_name = save_name+'.png'
plt.savefig(plot_name)
# plt.savefig(result_dir + plot_name
plt.show()  # Need to comment out plt.show to successfully save image
