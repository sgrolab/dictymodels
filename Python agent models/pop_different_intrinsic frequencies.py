# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 12:38:02 2019

@author: Chuqiao

Population of cells with different intrinsic frequencies

"""

import pandas as pd
import numpy as np
import random
import math
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import find_peaks
from time import perf_counter 

Nt_Sgro = 27
#%% Single cells oscillatory regime
from Sgro2015_agent_and_pop_FUN import Sgro2015_agent

e=0.02; tauA=0.09; tauR=tauA/e; g=0.5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
# Nt=27; # normalization factor of t
Na=3.5;  # normalization factor of A
A0=-1.5; R0=-0.5
Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)
Nt_Sgro = 27
dt=0.005 ; t_tot=20*Nt_Sgro; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
#constant_signal=1 
constant_signal=10000
stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
signal_trace=np.zeros(len(t))
signal_trace[stim_time_step:] = constant_signal


# test update function with one step caculations
# A,R,r=test_agent.update(1,dt)

# initializations
A_trace_orig=[A0]; R_trace_orig=[R0]; r_trace=[]

for i in range(len(t)-1):
    A_now=A_trace_orig[i]
    R_now=R_trace_orig[i]
    signal_now=signal_trace[i]
    
    A_next,R_next,r_now=Sgro_agent.update(dt,signal_now)
    A_trace_orig.append(A_next)
    R_trace_orig.append(R_next)
    r_trace.append(r_now)
    
# Traces
A_trace_offset=1.5
A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
A_trace_plot=(A_trace_orig+A_trace_offset)/Na;

## Check find_peaks
peaks, properties = find_peaks(A_trace_plot,  prominence=0.95)

#plt.plot(A_trace_plot)
#plt.plot(peaks, A_trace_plot[peaks], "x")
#Wdths = peak_widths(A_trace_plot, peaks, rel_height=0.95)
#plt.hlines(*Wdths[1:], color="C2")
#plt.show()
#Nt_Sgro = round(Wdths[0][0] * dt,2)
#print('Time normalization factor for Sgro 2015 is '+str(Nt_Sgro))
t_plot_Sgro = np.array(t)/(Nt_Sgro)
OscFreq = len(peaks)/t_plot_Sgro[-1]
label_font_size = 12
trace_width = 3

fig5 = plt.figure(figsize=(6, 8))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.3)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Sgro,signal_trace, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
ax1.set_title('cAMP stim:'+str(constant_signal))
# ax1.set_xlabel('Time',fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot_Sgro,A_trace_plot, color='g',linewidth=trace_width)
ax2.set_ylabel('Activator',fontsize=label_font_size)
ax2.yaxis.label.set_color('g')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax2.set_xlabel('Time',fontsize=label_font_size)
fig5.text(0.5, 0.02, 'Sgro single cell, e='+str(e)+',\n Osc frequency='+str(OscFreq), fontsize=label_font_size, ha='center')
plt.show()
#%%Population oscillations
from Sgro2015_agent_and_pop_FUN import Sgro2015_pop_mixed_cells
g=0.5; sigma = 0.15 # noise strength
N = 100 # number of cells in the population
rho = 10**(-3.5); j = 0.5

perc_WT_cells = 50 # percentage of normal cells
e_prime = 0.005
e=0.1*np.ones(N); e[int(perc_WT_cells/100*N):] = e_prime


SgroPopParam={'e':e,'g':g,'c0':1.2,'sigma':sigma,'N':N,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0, 'rho': rho,'j': j}
Na=3.5;  # normalization factor of A
A0=-1.5*np.ones(N) # + np.random.uniform(-sigma,sigma,N); ###########
R0=-0.5*np.ones(N) # + np.random.uniform(-sigma,sigma,N); 
cAMPext0 = 0
Sgro_pop_mixed_cells=Sgro2015_pop_mixed_cells(A0,R0,cAMPext0, SgroPopParam)

dt=0.005 ; t_tot=50*Nt_Sgro; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
alphafval =0
time_separation = 0

# initializations
A_trace_orig=np.zeros((N,len(t)))
R_trace_orig=np.zeros((N,len(t)))
cAMPext_trace = np.zeros((len(t),1))
A_trace_orig[:,0]=A0 
R_trace_orig[:,0]=R0
cAMPext_trace[0] = cAMPext0

for i in range(len(t)-1):
    A_now=A_trace_orig[:,i]
    R_now=R_trace_orig[:,i]
    cAMPext_now = cAMPext_trace[i]
    
    A_next,R_next, cAMPext_next = Sgro_pop_mixed_cells.update( dt, time_separation,alphafval)
    A_trace_orig[:,i] = A_next
    R_trace_orig[:,i] = R_next
    cAMPext_trace[i] = cAMPext_next

    
# Traces
A_trace_offset=1.5
A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
A_trace_plot=(A_trace_orig+A_trace_offset)/Na;
t_plot_Sgro = np.array(t)/Nt_Sgro

A_trace_mean_plot = np.mean(A_trace_plot,axis = 0)
A_trace_mean_plot_WT = np.mean(A_trace_plot[0:int(perc_WT_cells/100*N),:],axis = 0)
A_trace_mean_plot_low_e = np.mean(A_trace_plot[int(perc_WT_cells/100*N):,:],axis = 0)

cAMPext_trace_plot = np.array(cAMPext_trace)

#%  check simulation traces
label_font_size=18; trace_width=3; tick_font_size=18
fig3 = plt.figure(figsize=(12,12))
grid = plt.GridSpec(4, 1, wspace=0.2, hspace=0.45)

ax1= fig3.add_subplot(grid[0, 0])
ax1.plot(t_plot_Sgro,A_trace_plot[2:10,:].T, linewidth=1, label= r'activator, $cAMP_{cyto}$')
ax1.plot(t_plot_Sgro,A_trace_mean_plot, color='lightgreen', linewidth=trace_width)
# ax.set_ylim([-0.2,1.3])
#ax1.set_xlabel('Time')
#ax1.set_ylabel(r'Activator, $cAMP_{cyto}$')
ax1.set_title('Population mean')

ax2= fig3.add_subplot(grid[1, 0])
ax2.plot(t_plot_Sgro,A_trace_plot[2:10,:].T, linewidth=1, label= r'activator, $cAMP_{cyto}$')
ax2.plot(t_plot_Sgro,A_trace_mean_plot_WT, color='lightgreen', linewidth=trace_width)
# ax.set_ylim([-0.2,1.3])
ax2.set_title('WT mean')

ax3= fig3.add_subplot(grid[2, 0])
ax3.plot(t_plot_Sgro,A_trace_plot[85:95,:].T, linewidth=1, label= r'activator, $cAMP_{cyto}$')
ax3.plot(t_plot_Sgro,A_trace_mean_plot_low_e, color='lightgreen', linewidth=trace_width)
# ax.set_ylim([-0.2,1.3])
ax3.set_title('e='+str(e_prime)+ ' mean')

ax4= fig3.add_subplot(grid[3, 0])
ax4.plot(t_plot_Sgro,cAMPext_trace_plot, color='lightgreen', linewidth=trace_width)
# ax.set_ylim([-0.2,1.3])
ax4.set_title('extracellular cAMP')


fig3.text(0.5, 0.93, 'Sgro 2015 population oscillation,'+str(perc_WT_cells)+'% WT cells,\n'+r'$rho$= '+str(rho)+
          ', j= '+str(j)+',time separation= '+str(time_separation), fontsize=label_font_size, ha='center')
fig3.text(0.04, 0.5, r'$cAMP_{cyto}$', fontsize=label_font_size,va='center', rotation='vertical')
fig3.text(0.5, 0.04, 'Time, A.U.', fontsize=label_font_size, ha='center')
plt.show()

## Get the extracellular cAMP oscillation height
#later_portion = 0.2 # start count peaks after this X total simulation time
#cAMPext_trace_plot_later=cAMPext_trace_plot[math.floor(len(t)*later_portion):] # the later part of trace
#PkPos, PkProperties = find_peaks(cAMPext_trace_plot_later[:,0], prominence=(0.5,500))
#plt.figure()
#plt.plot(cAMPext_trace_plot_later)
#plt.plot(PkPos, cAMPext_trace_plot_later[PkPos], "x")

## Get the oscillation period
#later_portion = 0.2 # start count peaks after this X total simulation time
#A_trace_mean_plot_later=A_trace_mean_plot[math.floor(len(t)*later_portion):] # the later part of trace
#PkPos, PkProperties = find_peaks(A_trace_mean_plot_later, prominence=(0.5,1.5))
## Check find_peaks
#plt.plot(A_trace_mean_plot_later)
#plt.plot(PkPos, A_trace_mean_plot_later[PkPos], "x")
#
#Sgro_pop_osc_period = (1-later_portion)*t_tot / len(PkPos)
#print('group oscillation period for Sgro is '+str(Sgro_pop_osc_period))

