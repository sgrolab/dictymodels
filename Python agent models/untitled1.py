# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 11:15:45 2020

@author: ellin
"""


## Notebook Sgro with noise test

import numpy as np
import random
import math
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import chirp, find_peaks, peak_widths
import pandas as pd
import scipy.io

from Gregor2010_agent_and_pop_FUN import  Gregor2010_agent
from Sgro2015_agent_and_pop_FUN import Sgro2015_agent
from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_agent_3var
from MaedaLoomis2004_agent_and_pop_FUN import MaedaLoomis2004_agent
from Kamino2017_agent_and_pop_FUN import Kamino2017_agent 

from NB_SC_functions import * 
from NormParam import *
#%%
e=0.1; tauA=0.09; tauR=tauA/e; g=0.5
sigma = 0.15 # noise strength
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':sigma,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
dt=0.005 ; t_tot=30*Nt_Sgro; t=list(np.arange(0,t_tot,dt))
cAMP = 1
#%%
z0First_space_Sgro = np.array([1,2,4,8])
# fold change of the secondary spike compared to priming concentration # np.logspace(0.5, 0.5, num=1) 
FC_space_Sgro = np.logspace(0.1, 1.6, num=8)
prm_lims_Sgro = [0.4,2]
# Define the time of step stimulation is applied  
stim_time_step1 = int(round(1/3*t_tot/dt))
stim_time_step2=int(round(2/3*t_tot/dt))

num_of_runs = 3
# define which traces to plot, 1st column defines the index of priming concentration, 
# 2nd column defines the index of fold change
single_trace_to_plot = np.array([[0,1],[0,3],[0,5]])
PkPrm_noise_Sgro,  PkPrm_noise_Sgro_norm =SC_FCD_Noise(z0First_space_Sgro, FC_space_Sgro, 
                                                                      num_of_runs,cAMP, Nt_Sgro, dt, t,
                                                                      prm_lims_Sgro,stim_time_step1,stim_time_step2, single_trace_to_plot,
                                                                      Sgro2015_SC, SgroAgentParam, Nh_Sgro, Nh_Sgro_offset)

PkPrm_mean_noise_Sgro=np.mean(PkPrm_noise_Sgro,axis=2)   
PkPrm_se_noise_Sgro = np.std(PkPrm_noise_Sgro,axis=2)/math.sqrt(num_of_runs)

#%%
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Sgro)+1))
fig3 = plt.figure(figsize=(8, 7))
grid = plt.GridSpec(1, 1, wspace=0.1, hspace=0.35)
ax3= fig3.add_subplot(grid[0, 0])
for i in range(len(z0First_space_Sgro)):
    ax3.plot(FC_space_Sgro, PkPrm_mean_noise_Sgro[i,:],'o-', color = colors[i], lw = 3, ms = 10,  label='Priming cAMP='+str(z0First_space_Sgro[i]))
    ax3.errorbar(FC_space_Sgro, PkPrm_mean_noise_Sgro[i,:], yerr=PkPrm_se_noise_Sgro[i,:],
                 fmt = 'o', color=colors[i], ecolor= colors[i], elinewidth=3, capsize=10, capthick=3)
ax3.set_ylim([-0.1,1])
ax3.set_ylabel( '2nd spike prominence',fontsize=20)
ax3.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=20)
ax3.set_xscale('log')
ax3.tick_params(axis='both', which='major', labelsize=16)
ax3.set_title('Sgro 2015 (with noise)', fontdict={'fontsize': 20, 'fontweight': 'medium'})
leg = ax3.legend();
ax3.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': 15})