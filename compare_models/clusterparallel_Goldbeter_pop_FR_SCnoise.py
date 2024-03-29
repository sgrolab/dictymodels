# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 18:05:24 2019

@author: Chuqiao
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
#%% Golbeter 1987
from Goldbeter1987_agent_and_pop_FUN  import Goldbeter1987_pop_3var_SCnoise
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

sigma = 0 # noise strength
N = 100  # number of cells

#kc=5.4 # 3.58 #j
#h=5 #rho
ParamArrLen = 25
# Set up parameter sweep matrix for dilution rate and density
kc_arr = np.linspace(10**(-0.5), 10**(2), num=ParamArrLen)
one_over_h_arr = np.linspace(0.01, 10, num=ParamArrLen) # cell density array
h_arr = 1/one_over_h_arr 

# kc_arr = np.linspace(40.5,44.5,num = 4)
# one_over_h_arr = np.array([1])# cell density array
# h_arr = 1/one_over_h_arr 

Goldbeter3PopParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':0,'h':0, 'sigma':sigma, 'N':N}

#update function
dt=0.001; t_tot=30*Nt_Goldbeter; t=list(np.arange(0,t_tot,dt))
nSteps = len(t)
t_plot_Goldbeter = np.array(t); t_plot_Goldbeter = t_plot_Goldbeter/Nt_Goldbeter
signal_input = 0
pk_find_prm = 0.1
pk_find_thresh = 5/Nh_Goldbeter

def calc_updates_Goldbeter(kc_arr, h_arr, Goldbeter3PopParam, nSteps, index):
    ''' index: (kc_index,h_index) '''
    Goldbeter3PopParam['kc'] = kc_arr[index[0]]
    Goldbeter3PopParam['h'] =  h_arr[index[1]]
    

    p0=0.8; a0=3; b0=0.9; g0=0
    p_trace=np.zeros((N,nSteps));p_trace[:,0] = p0; 
    b_trace=np.zeros((N,nSteps));b_trace[:,0] = b0; 
    g_trace=np.zeros((nSteps,1));g_trace[0] = p0; 

    Goldbeter3_pop=Goldbeter1987_pop_3var_SCnoise(0,p0,b0,g0,Goldbeter3PopParam)
    for i in range(len(t)-1):
        p_next,b_next,g_next= Goldbeter3_pop.update(dt,a0,signal_input)
        
        p_trace[:,i+1]= p_next
        b_trace[:,i+1]= b_next
        g_trace[i+1] = g_next
        
   
    b_trace_mean = np.mean(b_trace,axis = 0)
    b_trace_mean_norm = b_trace_mean/Nh_Goldbeter
    later_portion = 0.2 # start count peaks after this X total simulation time
    b_trace_mean_norm_later=b_trace_mean_norm[math.floor(nSteps * later_portion):] # the later part of trace
    
    
    #p_trace_later=p_trace[math.floor(nSteps * later_portion):]
    #p_trace_later_norm = p_trace_later/np.amax(p_trace_later)
    #g_trace_later=g_trace[math.floor(nSteps * later_portion):]
    #g_trace_later_norm = g_trace_later/np.amax(g_trace_later)
    t_plot_Goldbeter_later = t_plot_Goldbeter[math.floor(nSteps * later_portion):]

    
    # pop_max = np.amax(b_trace_mean_norm); pop_min = np.amin(b_trace_mean_norm)
    # PkPos, PkProperties = find_peaks(b_trace_mean_norm_later, prominence = pk_find_prm*(pop_max-pop_min))

    PkPos, PkProperties = find_peaks(b_trace_mean_norm_later, prominence = pk_find_thresh)
#    # Check find_peaks
#    fig = plt.figure()
#    plt.plot(b_trace_later)
#    plt.plot(PkPos, b_trace_later[PkPos], "x")
    
    if len(PkPos) == 0:
        firing_rate = 0; height = 0
    else: 
        firing_rate = len(PkPos)/(t_tot/Nt_Goldbeter*(1-later_portion))
        height = np.mean(PkProperties["prominences"])
   
#    # check simulation traces
#    fig = plt.figure(figsize=(10,2.5)); grid = plt.GridSpec(1, 2,hspace= 0.3)
#    ax1= fig.add_subplot(grid[0, 0])
#    ax1.plot(t_plot_Goldbeter_later,b_trace_later, color = 'g', linewidth = 3, label = 'cAMPi')
#    ax1.plot(t_plot_Goldbeter_later,p_trace_later, color ='b',label = 'Proportion of active receptor')
#    ax1.plot(t_plot_Goldbeter_later,g_trace_later, color ='k', label = 'cAMPe')
#    # ax1.set_xlim([10,20])
#    ax1.set_xlabel('Time')
#    ax1.set_ylabel('b')
#    ax1.set_title('dilutio rate= '+ '{:#.3n}'.format(np.float64(kc_arr[index[0]])) +
#    ', density= '+'{:#.3n}'.format(np.float64(one_over_h_arr[index[1]])) + 
#        ', FR = '+'{:#.3n}'.format(np.float64(firing_rate)))
#    leg = ax1.legend()
#                
#    ax2= fig.add_subplot(grid[0, 1])
#    ax2.plot(t_plot_Goldbeter_later,b_trace_later_norm, color = 'g', linewidth = 3)
#    ax2.plot(t_plot_Goldbeter_later,p_trace_later, color ='b')
#    ax2.plot(t_plot_Goldbeter_later,g_trace_later_norm, color ='k')
#    ax2.set_xlim([10,20])
#    ax2.set_xlabel('Time')
#    ax2.set_ylabel('Normalized trace, A.U.')
#    ax2.set_title('Normalized traces, zoomed in')
#
#    plt.show()
    
    return index,firing_rate,height
    

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

#%% plot heat map


# assign results to pop_rate matrix
save_name = 'pop_FR_Goldbeter_200416_hnorm_dt'+str(dt)+'_noise'+str(sigma)+ \
'ParamLen'+str(ParamArrLen)+ '_pkfindthr'+str(pk_find_thresh)
pop_rate_Goldbeter = np.zeros([kc_arr.size, h_arr.size])
pop_height_Goldbeter = np.zeros([kc_arr.size, h_arr.size])
for i in results:
    idx = i[0]; firing_rate = i[1]; height = i[2]
    pop_rate_Goldbeter[idx] = firing_rate
    pop_height_Goldbeter[idx] = height
np.savez(save_name+'.npz', kc_arr = kc_arr, h_arr = h_arr,pop_rate_Goldbeter = pop_rate_Goldbeter, pop_height_Goldbeter=pop_height_Goldbeter)
    
# Plot heat map
title_font_size = 20 
label_font_size = 20
tick_font_size = 16

fig3 = plt.figure(figsize=(23, 4.5))
grid = plt.GridSpec(1,3, wspace=0.3, hspace=0.05)

ax1= fig3.add_subplot(grid[0,0])
heatmap = ax1.pcolor(kc_arr, one_over_h_arr, pop_rate_Goldbeter.transpose(), cmap='jet') # cmap='jet'
ax1.set_yscale('log')
# heatmap.set_clim(0,1.5)
cbar=fig3.colorbar(heatmap, ax=ax1);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}j$',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('dt='+str(dt)+', noise= '+str(sigma)+'Goldbeter1987 pop firing rate', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax2= fig3.add_subplot(grid[0,1])
heatmap = ax2.pcolor(kc_arr, one_over_h_arr, pop_height_Goldbeter.transpose(), cmap='jet') # cmap='jet'
ax2.set_yscale('log')
# heatmap.set_clim(0,2)
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}(\gamma)$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Goldbeter1987 pop firing height', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# save image and results to the current folder
plot_name = save_name+'.png'
plt.savefig(plot_name)

plt.show() # after saving image to successfully save the images




