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
## need to have the following to save images on the cluster
#if 'DISPLAY' not in os.environ:
#	matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Gregor2010_agent_and_pop_FUN import Gregor2010_pop
# Time scale normalization
Nt_Gregor = 6
Nt_Sgro = 27
Nt_Goldbeter = 7
Nt_Maeda = 3.57
Nt_Kamino = 6
#%% Gregor 2010

#set parameters
Amax=20;  Abas=0.4 # uM
w=2*math.pi/6 # min-1
Vc=1.1e-9 # ml
St=1.33 # cm2
Sc=1.3e-6 # cm2
K=0.0004 # uM, 400 pM
c_sec= 3.6 # min-1
c_excite=1.01 # min-1
Nc=100 # Num of cells
Vt = 1 #chamber size ml

ext_input = 0
time_separation = 0
eta=0.002 # noise stength

##create parameter arrays
#rho_arr = np.logspace(-3.5,1,num=3) # np.array([1e-2,1e-1,1,5,10]) #
#k_arr = np.linspace(1,50,num=3) # np.array([25]) #

k_arr=np.array([25]) 
rho_arr=np.array([1])

# Pack parameters together
GregorPopParam={'Amax':Amax,'Abas':Abas,'w':w,'Vc':Vc,'St':St,'Sc':Sc,'K':K,\
            'c_sec':c_sec,'c_excite':c_excite,'Nc':Nc}

dt=0.005; t_tot=30*Nt_Gregor; t=list(np.arange(0,t_tot,dt))
nSteps = len(t)

def calc_updates_Gregor(k_arr, rho_arr, GregorPopParam, nSteps, index):
    ''' index: (j_index,rho_index) '''
    k = k_arr[index[0]]
    rho = rho_arr[index[1]]
    
    campCyto0 = 7.5*np.ones(Nc)
    sinthetai0 = (campCyto0*2-Amax-Abas)/(-Amax+Abas)
    thetai0 = np.arcsin(sinthetai0)
    campExt0 = 0 # Vc*St/Sc*rho/K*c_sec*1/Nc*np.sum(campCyto0);

    Gregor_pop=Gregor2010_pop(campCyto0, thetai0, campExt0, GregorPopParam)
    
    gregor_thetai_trace=np.zeros((Nc,len(t))); gregor_thetai_trace[:,0] = thetai0
    gregor_campCyto_trace=np.zeros((Nc,len(t))); gregor_campCyto_trace[:,0] = campCyto0
    gregor_campExt_trace=np.zeros(len(t)); gregor_campExt_trace[0] = campExt0
    
    for i in range(nSteps-1):
        thetai_next, campCyto_next, campExt_next = Gregor_pop.update(dt,eta,rho,k,Vt,time_separation,ext_input)
        gregor_thetai_trace[:,i+1] = thetai_next
        gregor_campCyto_trace[:,i+1] = campCyto_next
        gregor_campExt_trace[i+1] = campExt_next
       
    #Traces
    gregor_campCyto_trace= np.array(gregor_campCyto_trace) 
    gregor_campCyto_trace_mean= np.mean(gregor_campCyto_trace,axis = 0)
    # gregor_campExt_trace = np.array(gregor_campExt_trace)

    later_portion = 0.2 # start count peaks after this X total simulation time
    gregor_campCyto_trace_mean_later=gregor_campCyto_trace_mean[math.floor(nSteps * later_portion):] # the later part of trace
    
    PkPos, PkProperties = find_peaks(gregor_campCyto_trace_mean_later, prominence=(0.5,30))

  
    if len(PkPos) == 0:
        firing_rate = 0; height = 0
    else: 
        firing_rate = len(PkPos)/(t_tot/Nt_Gregor*(1-later_portion))
        height = np.mean(PkProperties["prominences"])
    
    SC_traces = gregor_campCyto_trace;
    t_plot = np.array(t)/Nt_Gregor; cAMPi_mean = gregor_campCyto_trace_mean
    SC_traces_idx = [1,14,5,16,17,18] 
    # plot the traces
    fig = plt.figure(figsize=(6, 4))
    grid = plt.GridSpec(1, 1, wspace=0.5, hspace=0.3)
    ax1= fig.add_subplot(grid[0, 0])
    # plot defined single cell traces
    if SC_traces_idx != 0:
        for this_idx in SC_traces_idx:
            this_trace = SC_traces[this_idx,:]# /np.amax(SC_traces[this_idx,:])
            ax1.plot(t_plot,this_trace, color='b',alpha=0.6, linewidth=1.5)     
    # Plot population mean
    ax1.plot(t_plot, cAMPi_mean, color='g' ,linewidth=2, label = 'cAMPi population mean')
    ax1.set_title('Dilution rate= '+ '{:#.3n}'.format(np.float64(k)) +
        ', density= '+'{:#.3n}'.format(np.float64(rho)) +
            ',\n FR = '+'{:#.3n}'.format(np.float64(firing_rate)), fontdict={'fontsize':18})
    ax1.set_xlabel('Time, A.U.', fontsize=20); 
    ax1.set_ylabel('Raw traces',fontsize=20 )
    ax1.tick_params(grid_linewidth = 15, labelsize = 20)
    ax1.set_xlim([0, 30])
    leg = ax1.legend()
#    ax1.set_xlim([0,30]); ax1.set_ylim([-0.25,1.25])
    plt.show()
    
    return index,firing_rate, height

# To run the simulation get all indices for gamma and rho
indices = all_indices(k_arr, rho_arr)

# Create a partial function call to calc_updates
par_calc_updates_Sgro = functools.partial(calc_updates_Gregor, k_arr, rho_arr, GregorPopParam, nSteps)

n_workers = get_n_workers()

# Start up the worker pool and return the results as a list of
# [((gamma_index, rho_index), len(PkPos)),...]
tic = perf_counter()
results = list(map(par_calc_updates_Sgro,indices))

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
pop_rate_Gregor = np.zeros([k_arr.size, rho_arr.size]) # population firing rate
pop_height_Gregor = np.zeros([k_arr.size, rho_arr.size]) 

for i in results:
    idx = i[0]; firing_rate = i[1]; height = i[2]
    pop_rate_Gregor[idx] = firing_rate
    pop_height_Gregor[idx] = height
np.savez('pop_fire_rate_Gregor_OUT_191026.npz', k_arr = k_arr, rho_arr = rho_arr,pop_rate_Gregor = pop_rate_Gregor, pop_height_Gregor=pop_height_Gregor)
    
# Plot heat map
title_font_size = 16
label_font_size = 16
tick_font_size = 16

fig3 = plt.figure(figsize=(23, 4.5))
grid = plt.GridSpec(1,3, wspace=0.3, hspace=0.05)

ax1= fig3.add_subplot(grid[0,0])
heatmap = ax1.pcolor(k_arr, rho_arr, pop_rate_Gregor.transpose(), cmap='jet') # cmap='jet'
ax1.set_yscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax1);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}j$',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Gregor 2010 pop firing rate', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax2= fig3.add_subplot(grid[0,1])
heatmap = ax2.pcolor(k_arr, rho_arr, pop_height_Gregor.transpose(), cmap='jet') # cmap='jet'
ax1.set_yscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}(\gamma)$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Gregor 2010 pop firing height, TS='+str(time_separation)+', noise = '+str(eta),
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# save image and results to the current folder
plot_name = 'pop_fire_rate_Gregor_191026.png'
plt.savefig(plot_name)
plt.show()

