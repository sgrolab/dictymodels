# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 15:40:02 2020

@author: cqhuyan
figure4 concept demonstration figures

"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import chirp, find_peaks, peak_widths
import pandas as pd
import scipy.io

# Normalization parameters
from NormParam import *

# set up new default font
import matplotlib
font = {'family' : 'Arial'}
matplotlib.rc('font', **font)

from celluloid import Camera
#%%
mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']  
abcd_font_size = 28
label_font_size=24
title_font_size = 26
sublabel_font_size = 22
trace_width=3
tick_font_size=20

#%% Sgro 2015
from Sgro2015_agent_and_pop_FUN import Sgro2015_agent

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5; c0 = 1.2
a= 0.058; Kd = 1e-5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'flux_thrs':0}

A0=-1.5; R0=-0.5
dt=0.005 ; t_tot=6*Nt_Sgro; t=list(np.arange(0,t_tot,dt))
signal_1=1
stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
signal_trace_1=np.zeros(len(t))
signal_trace_1[stim_time_step:] = signal_1

A_trace_orig_1=[A0]; R_trace_orig_1=[R0]; r_trace_1=[]

# adaptive spike
Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)
for i in range(len(t)-1):
    A_now=A_trace_orig_1[i]
    R_now=R_trace_orig_1[i]
    signal_now_1=signal_trace_1[i]
    
    A_next,R_next,r_now=Sgro_agent.update(dt,signal_now_1)
    A_trace_orig_1.append(A_next)
    R_trace_orig_1.append(R_next)
    r_trace_1.append(r_now)    

A_trace_offset=1.5
A_trace_orig_1 = np.array(A_trace_orig_1) # vectorize A_trace_orig
R_trace_orig_1 = np.array(R_trace_orig_1)
A_trace_plot_1=(A_trace_orig_1+A_trace_offset)/Nh_Sgro;
# nullclines
A_null = np.linspace(-2.5,2.5,num=200)
dAdt_null_1=A_null-1/3*A_null**3+a*np.log(1+signal_1/Kd)
dAdt_null_no_stim=A_null-1/3*(A_null**3)
dRdt_null=1/g*(A_null+c0)


#%% Kamino excitability
from Kamino2017_agent_and_pop_FUN import Kamino2017_agent

tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
AgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}

dt=0.001 ; t_tot=6*Nt_Kamino; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
#constant_signal=1 
signal_1=1
stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
signal_trace_1=np.zeros(len(t))
signal_trace_1[stim_time_step:] = signal_1
# initializations
x0=0.01; y0=0.05; z0=0.005
x_trace_1=[x0]; y_trace_1=[y0]

# adaptive spike
Kamino_agent=Kamino2017_agent([x0,y0,z0],AgentParam)
for i in range(len(t)-1):
    x_now=x_trace_1[i]
    y_now=y_trace_1[i]
    x_next,y_next,z_next= Kamino_agent.update(dt, signal_trace_1[i])
    x_trace_1.append(x_next)
    y_trace_1.append(y_next)

x_trace_1 = np.array(x_trace_1) 
y_trace_1 = np.array(y_trace_1)
y_trace_1_hnorm = (y_trace_1-Nh_Kamino_offset)/Nh_Kamino

t_plot_Kamino = np.array(t)/(Nt_Kamino)

# nullclines
x_null_short = np.linspace(-0.5,2,num=100)
x_null_long = np.linspace(-1000,11000,num=100)
dxdt_null_no_stim = 0 + delta
dxdt_null_1 = signal_1 + delta
dydt_null_no_stim_short = (0+delta)**n/((0+delta)**n+(np.power(K*x_null_short,n)))
dydt_null_no_stim_long = (0+delta)**n/((0+delta)**n+(np.power(K*x_null_long,n)))
dydt_null_1 = (signal_1+delta)**n/((signal_1+delta)**n+(np.power(K*x_null_short,n)))




#%% plot excitability
abcd_font_size = 28
label_font_size=24
title_font_size = 26
sublabel_font_size = 22
trace_width=2
tick_font_size=20

fig5 = plt.figure(figsize=(16,5))
grid = plt.GridSpec(3,3, wspace=0.25, hspace=0.1)

ax1 = fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Kamino, signal_trace_1, 'k',linewidth=trace_width)
ax1.set_ylabel( 'Signal' ,fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_xlim([0,6]);
ax1.set_xticks([]); ax1.set_yticks([]);

ax0 = fig5.add_subplot(grid[1:, 0])
ax0.plot(t_plot_Kamino, y_trace_1_hnorm, 'k',linewidth=trace_width)
# ax0.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
# ax0.set_title(r'$cAMP_{e}$ input: '+str(signal_1),fontsize=label_font_size)
# ax0.tick_params(axis='both', which='both', bottom=False, top=False,labelbottom=False,labelsize=tick_font_size)
ax0.set_xticks([]); ax0.set_yticks([]);
ax0.set_xlabel( 'Time' ,fontsize=label_font_size)
ax0.set_ylabel( 'Output' ,fontsize=label_font_size)
ax0.set_xlim([0,6]); ax0.set_ylim([-0.15,1.05])

# Kamino
ax2= fig5.add_subplot(grid[:, 1])
ax2.axvline(x=dxdt_null_no_stim, ls='-', linewidth=trace_width, color='blue', linestyle = '--')
ax2.axvline(x=dxdt_null_1, ls='-', linewidth=trace_width, color='blue')
ax2.plot(x_null_short,dydt_null_no_stim_short,'firebrick',linewidth=trace_width)
ax2.plot(x_null_short, dydt_null_1,'firebrick',linestyle = '--',linewidth=trace_width)

ax2.plot( x_trace_1, y_trace_1,'k',linewidth=trace_width)
ax2.plot( x_trace_1[0], y_trace_1[0],'o',markersize =7,color='k')
ax2.plot( x_trace_1[-1], y_trace_1[-1],'o',markersize =7,color='k',fillstyle='none')
ax2.set_xlim([-0.2,1.2])
ax2.set_ylim([-0.05,1.05])
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)

# vector field for cAMPe input = 1
x_arr = np.arange(-1,3,0.1)
y_arr = np.arange(-0.5,1.5,0.1)
x_mesh, y_mesh = np.meshgrid(x_arr , y_arr )
makelong = 1
dx = 1/tau*(signal_1+delta-x_mesh); dx = makelong*dx
dy = (signal_1+delta)**n/((signal_1+delta)**n+(np.power(K*x_mesh,n))) - y_mesh;  dy = makelong*dy
# ax2.streamplot(x_mesh,y_mesh,dx, dy,color='forestgreen')
q = ax2.quiver(x_mesh,y_mesh,dx, dy,scale=3, scale_units='inches',width=0.005,color='grey')
ax2.set_xticks([]); ax2.set_yticks([]);
ax2.set_xlabel( 'I' ,fontsize=label_font_size)
ax2.set_ylabel( 'A' ,fontsize=label_font_size)

#Sgro
ax3= fig5.add_subplot(grid[:, 2])
ax3.plot(A_null,dRdt_null, color='firebrick', linewidth=trace_width)
ax3.plot(A_null,dAdt_null_no_stim,'blue',linewidth=trace_width)
ax3.plot(A_null, dAdt_null_1,'blue',linestyle = '--', linewidth=trace_width)
ax3.plot( A_trace_orig_1, R_trace_orig_1,'k',linewidth=trace_width)
ax3.plot(A_trace_orig_1[0],R_trace_orig_1[0],'o',markersize = 7,color='k')
ax3.plot( A_trace_orig_1[-1],R_trace_orig_1[-1], 'o',markersize =7,color='k',fillstyle='none')
ax3.set_ylim([-1,2]);
ax3.set_xlim([-2.5,2.5])

ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)

A_arr = np.arange(-2.5,3,0.6)
R_arr = np.arange(-1.2,3,0.3)
makelong = 1
A_mesh, R_mesh = np.meshgrid(A_arr , R_arr )
dA = A_mesh-(np.power(A_mesh,3))/3-R_mesh+a*np.log(1+signal_1/Kd); dA=makelong*dA
dR = e*(A_mesh-g* R_mesh+c0); dR=makelong*dR
ax3.quiver(A_mesh,R_mesh,dA, dR,scale=9, scale_units='inches',width=0.005,color='grey')
ax3.set_xticks([]); ax3.set_yticks([]);
ax3.set_xlabel( 'A' ,fontsize=label_font_size)
ax3.set_ylabel( 'I' ,fontsize=label_font_size)


# ax1= fig5.add_subplot(grid[1:, 0])
# ax1.axvline(x=dxdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
# ax1.axvline(x=dxdt_null_1, ls='-', linewidth=trace_width, color='dimgrey')
# ax1.plot(x_null_short,dydt_null_no_stim_short,'lightblue',linewidth=trace_width)
# ax1.plot(x_null_short, dydt_null_1,'deepskyblue',linewidth=trace_width)

# ax1.plot( x_trace_1, y_trace_1,'green',linewidth=trace_width)
# ax1.plot( x_trace_1[0], y_trace_1[0],'*',markersize = 12,color='springgreen')
# ax1.set_xlim([-0.5,2])
# ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)



#%% FCD demonstration
# Kamino 
from Kamino2017_agent_and_pop_FUN import Kamino2017_agent

tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
AgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}

dt=0.001 ; t_tot=9*Nt_Kamino; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
#constant_signal=1 
signal_1=1; time_1=0# at this time step input 1 is applied
signal_2 = 3; time_2= int(round(3/9*t_tot/dt))
signal_3 = 9; time_3= int(round(6/9*t_tot/dt))
signal_trace=np.zeros(len(t))
signal_trace[time_1:] = signal_1
signal_trace[time_2:] = signal_2
signal_trace[time_3:] = signal_3

# initializations
x0=0.01; y0=0.05; z0=0.005
x_trace=[x0]; y_trace=[y0]

# adaptive spike
Kamino_agent=Kamino2017_agent([x0,y0,z0],AgentParam)
for i in range(len(t)-1):
    x_now=x_trace[i]
    y_now=y_trace[i]
    x_next,y_next,z_next= Kamino_agent.update(dt, signal_trace[i])
    x_trace.append(x_next)
    y_trace.append(y_next)

x_trace = np.array(x_trace) 
y_trace = np.array(y_trace)
y_trace_hnorm = (y_trace-Nh_Kamino_offset)/Nh_Kamino

t_plot_Kamino = np.array(t)/(Nt_Kamino)

#plot demonstration
abcd_font_size = 28
label_font_size=24
title_font_size = 26
sublabel_font_size = 22
trace_width=2
tick_font_size=20

fig5 = plt.figure(figsize=(5.5,5))
grid = plt.GridSpec(3,1, wspace=0.25, hspace=0.1)

ax1 = fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Kamino, signal_trace, 'k',linewidth=trace_width)
ax1.set_ylabel( 'Signal' ,fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_xlim([2,9]);
ax1.set_xticks([]); ax1.set_yticks((1,3,9));
ax1.tick_params(axis='both', which='both', labelsize=tick_font_size)

ax0 = fig5.add_subplot(grid[1:, 0])
ax0.plot(t_plot_Kamino, y_trace_hnorm, 'k',linewidth=trace_width)
# ax0.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
# ax0.set_title(r'$cAMP_{e}$ input: '+str(signal_1),fontsize=label_font_size)
# ax0.tick_params(axis='both', which='both', bottom=False, top=False,labelbottom=False,labelsize=tick_font_size)
ax0.set_xticks([]); ax0.set_yticks([]);
ax0.set_xlabel( 'Time' ,fontsize=label_font_size)
ax0.set_ylabel( 'Output' ,fontsize=label_font_size)
ax0.set_xlim([2,9]); ax0.set_ylim([-0.05,0.4])

