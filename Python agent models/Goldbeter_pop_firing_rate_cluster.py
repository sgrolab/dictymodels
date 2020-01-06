# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 18:05:24 2019

@author: ellin
"""
import numpy as np
from time import perf_counter 
import functools
import multiprocessing as mp
from scipy.signal import find_peaks
import math

from Population_parallel_tools import all_indices, get_n_workers
import matplotlib.pyplot as plt

# needs to be determined
Nt_Gregor = 6
Nt_Sgro = 27
Nt_Goldbeter = 7
Nt_Maeda = 3.57
Nt_Kamino = 6
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

kc_arr = np.linspace(0, 10, num=2)
h_arr = np.linspace(0, 10, num=2)

Goldbeter3PopParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':0,'h':0}

#update function
dt=0.001; t_tot=100*Nt_Goldbeter; t=list(np.arange(0,t_tot,dt))
nSteps = len(t)
t_plot_Goldbeter = np.array(t); t_plot_Goldbeter = t_plot_Goldbeter/Nt_Goldbeter
signal_input = 0

def calc_updates_Goldbeter(kc_arr, h_arr, Goldbeter3PopParam, nSteps, index):
    ''' index: (kc_index,h_index) '''
    Goldbeter3PopParam['kc'] = 11 #kc_arr[index[0]]
    Goldbeter3PopParam['h'] =  1#h_arr[index[1]]
    

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
    later_portion = 0 # start count peaks after this X total simulation time
    b_trace_later=b_trace[math.floor(nSteps * later_portion):] # the later part of trace
##     check simulation traces
#    fig = plt.figure()
#    plt.plot(t_plot_Goldbeter,b_trace)
#    plt.xlabel('Time')
#    plt.ylabel('b')
#    plt.title('kc= '+str(Goldbeter3PopParam['kc'])+' h= '+str(Goldbeter3PopParam['h']))
#    plt.show()

    PkPos, PkProperties = find_peaks(b_trace_later, prominence=(5,1000))
#    # Check find_peaks
#    fig = plt.figure()
#    plt.plot(b_trace_later)
#    plt.plot(PkPos, b_trace_later[PkPos], "x")
#    plt.title('kc= '+str(Goldbeter3PopParam['kc'])+' h= '+str(Goldbeter3PopParam['h']))
#    plt.show()
    
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
results = list(map(par_calc_updates_Goldbeter,indices))

## Run serially or in parallel if possible
#if n_workers == 1:
#    results = list(map(par_calc_updates_Goldbeter,indices))
#else:
#    with mp.Pool(get_n_workers()) as pool:
#        results = pool.map(par_calc_updates_Goldbeter,indices)
#
#toc = perf_counter() 
#print('time passed with %s workers: %s' % (get_n_workers(),toc-tic))


#%%
# assign results to pop_rate matrix
pop_rate_Goldbeter = np.zeros([kc_arr.size, h_arr.size])
for i in results:
    idx = i[0]; firing_rate = i[1]
    pop_rate_Goldbeter[idx] = firing_rate
    
# Plot heat map
title_font_size = 20 
label_font_size = 20
tick_font_size = 16

fig3 = plt.figure(figsize=(23, 4.5))
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

# save image and results to the current folder
# results_dir = 'C:/Users/ellin/Documents/GitHub/dictymodels/Python agent models/cluster parallel/'
plot_name = 'C:/Users/ellin/Documents/GitHub/dictymodels/Python agent models/pop_fire_rate_Goldbeter_test.png'
plt.savefig(plot_name)
# plt.savefig(result_dir + plot_name)

np.savez('pop_fire_rate_Goldbeter_OUT.npz', kc_arr = kc_arr, h_arr = h_arr,pop_rate_Goldbeter = pop_rate_Goldbeter)
## Load saved arrays
#data = np.load('pop_fire_rate_Goldbeter_OUT.npz')