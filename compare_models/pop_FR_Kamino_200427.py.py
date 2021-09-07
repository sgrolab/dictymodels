# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 22:29:46 2020

@author: cqhuyan

FCD of Kamino: is the increase in pop firing rate at higher dilution 
rates is an essential property of the FCD network?
"""

import os
import numpy as np
from time import perf_counter 
from scipy.signal import find_peaks
import math
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd

from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_pop_3var
from MaedaLoomis2004_agent_and_pop_FUN import MaedaLoomis2004_pop
from Gregor2010_agent_and_pop_FUN import Gregor2010_pop
from Sgro2015_agent_and_pop_FUN import Sgro2015_pop
from Kamino2017_agent_and_pop_FUN import Kamino2017_pop
from NB_pop_functions import * 
# Normalization parameters
from NormParam import *

#%%
tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3 # cAMPe dilution/degradation rate 
rho= 1 # population density
KaminoPopParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
dt=0.001; t_tot = 30; t=list(np.arange(0,t_tot*Nt_Kamino,dt))
nSteps = len(t)
#%%
# cAMPe influx level. Default is 0. 
cAMPext_influx_trace=np.zeros(len(t))

gamma_arr = np.array([3,5,7,10])# ([25]) # 
rho_arr = np.array([10]) # np.array([20,15,10]) # ([10])
if_plt_traces = True # whether to plot individuial traces

for gamma_temp in gamma_arr:
    KaminoPopParam['gamma'] = gamma_temp
    for rho_temp in rho_arr:
        KaminoPopParam['rho'] = rho_temp
        [t_plot_Kamino, y_trace, x_trace, z_trace] = Kamino_pop(KaminoPopParam,dt,t,cAMPext_influx_trace)
        t_plot_Kamino = t_plot_Kamino/(Nt_Kamino)
        if if_plt_traces:
            later_portion = 0.2 # start count peak after this X total simulation time
            y_trace_later=y_trace[math.floor(nSteps * later_portion):] # the later part of trace
            # y_trace_later_norm = y_trace_later/np.amax(y_trace_later)
            x_trace_later=x_trace[math.floor(nSteps * later_portion):] 
            # x_trace_later_norm = x_trace_later/np.amax(x_trace_later)
            z_trace_later=z_trace[math.floor(nSteps * later_portion):] 
            # z_trace_later_norm = z_trace_later/np.amax(z_trace_later)
            t_plot_Kamino_later = t_plot_Kamino[math.floor(nSteps * later_portion):]

            PkPos, PkProperties = find_peaks(y_trace_later, prominence=(0.02,100))
            if len(PkPos) == 0:
                firing_rate = 0; height = 0
            else: 
                firing_rate = len(PkPos)/(t_tot*(1-later_portion))
                height = np.mean(PkProperties["prominences"])
            l_base = PkProperties["left_bases"] [1]  
            pkpos = PkPos[1]
            FC = y_trace_later[pkpos]/y_trace_later[l_base]
            
            # check simulation traces
            fig = plt.figure(figsize=(6,5)); grid = plt.GridSpec(1, 1,hspace= 0.3,wspace = 0.5)
            ax1= fig.add_subplot(grid[0, 0])
            ax1.plot(t_plot_Kamino_later,y_trace_later, color = 'g', linewidth = 1, label = 'cAMPi')
            
            ax1.plot(t_plot_Kamino_later[l_base],y_trace_later[l_base], color = 'g', marker = '*')
            ax1.plot(t_plot_Kamino_later[pkpos],y_trace_later[pkpos], color = 'g', marker = 'x')
            
            ax1.plot(t_plot_Kamino_later,x_trace_later, color ='b',label = 'inhibitor')
            ax1.plot(t_plot_Kamino_later,z_trace_later, color ='k', label = 'cAMPe')
            
            # ax1.set_xlim([10,25])
            # ax1.set_ylim([0,3.5])
            ax1.set_xlabel('Time, A.U.', fontsize=20); 
            ax1.set_ylabel('Traces',fontsize=20 )
            ax1.tick_params(grid_linewidth = 15, labelsize = 20)
            ax1.set_title('Dilution rate= '+ '{:#.3n}'.format(np.float64(gamma_temp)) +
            ', density= '+'{:#.3n}'.format(np.float64(rho_temp)) + 
                '\nFR = '+'{:#.3n}'.format(np.float64(firing_rate))+
                ', FC = '+'{:#.3n}'.format(np.float64(FC)),
                         fontdict={'fontsize':18})
            leg = ax1.legend()
            # ax2= fig.add_subplot(grid[0, 1])
            # ax2.plot(t_plot_Kamino_later,x_trace_later_norm, color = 'g', linewidth = 3)
            # ax2.plot(t_plot_Kamino_later,y_trace_later, color ='b')
            # ax2.plot(t_plot_Kamino_later,z_trace_later_norm, color ='k')
            # ax2.set_xlim([30,50])
            # ax2.set_xlabel('Time')
            # ax2.set_ylabel('Normalized traces, A.U.')
            # ax2.set_title('Normalized traces, zoomed in')

            plt.show()