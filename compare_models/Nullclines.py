# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 11:02:34 2019

@author: Chuqiao 

Nullcline and phase portrait of Goldbeter 1987, Sgro 2015, and Kamino 2017

"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import chirp, find_peaks, peak_widths
import pandas as pd
import scipy.io
import sys

# Normalization parameters
sys.path.append("//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models")
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

# experimental ramp
my_dir = '//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/exp_data/'
Sgro2015Figure3excel = pd.read_excel(my_dir+'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure3')

#%% Sgro 2015 nullclines, adaptive spike and oscillations
from Sgro2015_agent_and_pop_FUN import Sgro2015_agent

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5; c0 = 1.2
a= 0.058; Kd = 1e-5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'flux_thrs':0}

A0=-1.5; R0=-0.5
dt=0.005 ; t_tot=6*Nt_Sgro; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
#constant_signal=1 
signal_1=1
signal_10k = 10000
stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
signal_trace_1=np.zeros(len(t))
signal_trace_1[stim_time_step:] = signal_1
signal_trace_10k=np.zeros(len(t))
signal_trace_10k[stim_time_step:] = signal_10k

# initializations
A_trace_orig_1=[A0]; R_trace_orig_1=[R0]; r_trace_1=[]
A_trace_orig_10k=[A0]; R_trace_orig_10k=[R0]; r_trace_10k=[]

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

# sustained oscillations
Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)
for i in range(len(t)-1):
    A_now=A_trace_orig_10k[i]
    R_now=R_trace_orig_10k[i]
    signal_now_10k=signal_trace_10k[i]
    
    A_next,R_next,r_now=Sgro_agent.update(dt,signal_now_10k)
    A_trace_orig_10k.append(A_next)
    R_trace_orig_10k.append(R_next)
    r_trace_10k.append(r_now)    

A_trace_offset=1.5
A_trace_orig_10k = np.array(A_trace_orig_10k) # vectorize A_trace_orig
R_trace_orig_10k = np.array(R_trace_orig_10k)
A_trace_plot_10k=(A_trace_orig_10k+A_trace_offset)/Nh_Sgro;

t_plot_Sgro = np.array(t)/(Nt_Sgro)
# nullclines
A_null = np.linspace(-2.5,2.5,num=200)
dAdt_null_no_stim=A_null-1/3*(A_null**3)
dAdt_null_1=A_null-1/3*A_null**3+a*np.log(1+signal_1/Kd)
dAdt_null_10k=A_null-1/3*A_null**3+a*np.log(1+signal_10k/Kd);
dRdt_null=1/g*(A_null+c0)

#%%
fig5 = plt.figure(figsize=(12,9))
grid = plt.GridSpec(3, 2, wspace=0.3, hspace=0.8)

ax0 = fig5.add_subplot(grid[0, 0])
ax0.plot(t_plot_Sgro, A_trace_plot_1, 'green',linewidth=trace_width)
ax0.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
ax0.set_title(r'$cAMP_{e}$ input: '+str(signal_1),fontsize=label_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_xlim([0,6]); ax0.set_ylim([-0.4,1.2])
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax1= fig5.add_subplot(grid[1:, 0])
ax1.plot(A_null,dRdt_null, 'darkgrey',linewidth=trace_width)
ax1.plot(A_null,dAdt_null_no_stim,'lightblue',linewidth=trace_width)
ax1.plot(A_null, dAdt_null_1,'deepskyblue',linewidth=trace_width)
ax1.plot( A_trace_orig_1, R_trace_orig_1,'green',linewidth=trace_width)
ax1.plot(A_trace_orig_1[0],R_trace_orig_1[0],'*',markersize = 12,color='springgreen')
ax1.set_ylim([-1,2.5])
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
#ax1.set_aspect('equal')

#cAMP_ext = 1
#A_arr = np.arange(-2.5,2.5,0.2)
#R_arr = np.arange(-1,2.5,0.2)
#A_mesh, R_mesh = np.meshgrid(A_arr , R_arr )
#dA = A_mesh-(np.power(A_mesh,3))/3-R_mesh+a*np.log(1+cAMP_ext/Kd)
#dR = e*(A_mesh-g* R_mesh+c0)
#ax1.quiver(A_mesh, R_mesh,dA,dR, linewidths=0.1, edgecolors='k')

ax00 = fig5.add_subplot(grid[0, 1])
ax00.plot(t_plot_Sgro, A_trace_plot_10k, 'green',linewidth=trace_width)
ax00.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
ax00.set_title(r'$cAMP_{e}$ input: '+str(signal_10k),fontsize=label_font_size)
ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax00.set_xlim([0,6]); ax00.set_ylim([-0.4,1.2])
ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax2= fig5.add_subplot(grid[1:, 1])
ax2.plot(A_null,dRdt_null, 'darkgrey',linewidth=trace_width)
ax2.plot(A_null,dAdt_null_no_stim,'lightblue',linewidth=trace_width)
ax2.plot(A_null, dAdt_null_10k,'deepskyblue',linewidth=trace_width)
ax2.plot( A_trace_orig_10k, R_trace_orig_10k,'green',linewidth=trace_width)
ax2.plot(A_trace_orig_10k[0],R_trace_orig_10k[0],'*',markersize = 12,color='springgreen')
ax2.set_ylim([-1,2.5])
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
#ax2.set_xlabel( 'Activator' ,fontsize=label_font_size)
#ax2.set_ylabel( 'Inhibitor' ,fontsize=label_font_size)
#ax2.set_aspect('equal')

#cAMP_ext = 10000
#A_arr = np.arange(-2.5,2.5,0.2)
#R_arr = np.arange(-1,2.5,0.2)
#A_mesh, R_mesh = np.meshgrid(A_arr , R_arr )
#dA = A_mesh-(np.power(A_mesh,3))/3-R_mesh+a*np.log(1+cAMP_ext/Kd)
#dR = e*(A_mesh-g* R_mesh+c0)
#ax2.quiver(A_mesh, R_mesh,dA,dR, linewidths=0.1, edgecolors='k')
fig5.text(0.04, 0.95,'B',ha='center', va='center', fontsize=abcd_font_size)

fig5.text(0.04, 0.8, r'$cAMP_{i}$'+'\n(Activator)', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)

fig5.text(0.5, 0.04,r'$cAMP_{i}$'+'(Activator)',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.04, 0.35, 'Inhibitor', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

fig5.text(0.5, 0.95,'Sgro 2015',ha='center', va='center',color=mycolors[5],fontsize=title_font_size)
plt.show()

# #%% save output in npz files
# np.savez('figS1_nullcline_200422.npz',
#          t_plot_Sgro = t_plot_Sgro,
#          signal_1 = signal_1, signal_10k = signal_10k,
#          A_trace_plot_1 = A_trace_plot_1,
#          A_trace_plot_10k = A_trace_plot_10k,
#          A_null = A_null, dRdt_null = dRdt_null,
#          dAdt_null_no_stim = dAdt_null_no_stim, 
#          dAdt_null_1 = dAdt_null_1, dAdt_null_10k = dAdt_null_10k,
#          A_trace_orig_1 = A_trace_orig_1 , R_trace_orig_1 = R_trace_orig_1,
#          A_trace_orig_10k = A_trace_orig_10k, R_trace_orig_10k = R_trace_orig_10k)

# #%% load data
# from iocustom import import_npz
# import_npz('figS1_nullcline_200422.npz',globals())
#%% nullclines make video
fig5 = plt.figure(figsize=(12,9))
grid = plt.GridSpec(3, 2, wspace=0.3, hspace=0.8)
fig5.text(0.04, 0.95,'B',ha='center', va='center', fontsize=abcd_font_size)

fig5.text(0.04, 0.8, r'$cAMP_{i}$'+'\n(Activator)', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)

fig5.text(0.5, 0.04,r'$cAMP_{i}$'+'(Activator)',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.04, 0.35, 'Inhibitor', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

fig5.text(0.5, 0.95,'Sgro 2015',ha='center', va='center',color=mycolors[5],fontsize=title_font_size)
ax0 = fig5.add_subplot(grid[0, 0])
ax1= fig5.add_subplot(grid[1:, 0])
ax00 = fig5.add_subplot(grid[0, 1])
ax2= fig5.add_subplot(grid[1:, 1])

camera = Camera(fig5)
down_sample_frame =  350

for huhu in range(int(len(t_plot_Sgro)/down_sample_frame)): 
    i = huhu*down_sample_frame
    ax0.plot(t_plot_Sgro[0:i], A_trace_plot_1[0:i], 'green',linewidth=trace_width)
    ax0.plot(t_plot_Sgro[i],A_trace_plot_1[i],'o',markersize = 8,color='springgreen')
    ax0.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
    ax0.set_title(r'$cAMP_{e}$ input: '+str(signal_1),fontsize=label_font_size)
    ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax0.set_xlim([0,6]); ax0.set_ylim([-0.4,1.2])
    ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    ax1.plot(A_null,dRdt_null, 'darkgrey',linewidth=trace_width)
    if i < stim_time_step:
        ax1.plot(A_null,dAdt_null_no_stim,'lightblue',linewidth=trace_width)
    else:
        ax1.plot(A_null, dAdt_null_1,'deepskyblue',linewidth=trace_width)
    ax1.plot( A_trace_orig_1[0:i], R_trace_orig_1[0:i],'green',linewidth=trace_width)
    ax1.plot(A_trace_orig_1[i],R_trace_orig_1[i],'o',markersize = 8,color='springgreen')
    ax1.set_ylim([-1,2.5])
    ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
    #ax1.set_aspect('equal')
    
    ax00.plot(t_plot_Sgro[:i], A_trace_plot_10k[:i], 'green',linewidth=trace_width)
    ax00.plot(t_plot_Sgro[i],A_trace_plot_10k[i],'o',markersize = 8,color='springgreen')
    ax00.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
    ax00.set_title(r'$cAMP_{e}$ input: '+str(signal_10k),fontsize=label_font_size)
    ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax00.set_xlim([0,6]); ax00.set_ylim([-0.4,1.2])
    ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    ax2.plot(A_null,dRdt_null, 'darkgrey',linewidth=trace_width)
    if i < stim_time_step:
        ax2.plot(A_null,dAdt_null_no_stim,'lightblue',linewidth=trace_width)
    else:
        ax2.plot(A_null, dAdt_null_10k,'deepskyblue',linewidth=trace_width)
    ax2.plot( A_trace_orig_10k[0:i], R_trace_orig_10k[0:i],'green',linewidth=trace_width)
    ax2.plot(A_trace_orig_10k[i],R_trace_orig_10k[i],'o',markersize = 8,color='springgreen')
    ax2.set_ylim([-1,2.5])
    ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    camera.snap()
    # plt.pause(0.05)
    # plt.clf()       
ani = camera.animate()

#%%
OutPath = 'C:\\Users\\ellin\\Dropbox\\Prospectus\\Presentation\\nullcline_videos\\'
ani_name = 'Sgro_bif'
ani.save(OutPath + ani_name + '.mp4') 
  
#%% phase portrait different excitability
from Sgro2015_agent_and_pop_FUN import Sgro2015_agent

e_small=0.1; e_big = 0.6
tauA=0.09; tauR=tauA/e; g=0.5; c0 = 1.2
a= 0.058; Kd = 1e-5
# Create agent with small or big e
SgroAgentParam={'e':e_small,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'flux_thrs':0}
Sgro_agent_smalle=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)
SgroAgentParam={'e':e_big,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'flux_thrs':0}
Sgro_agent_bige=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)

A0=-1.5; R0=-0.5
dt=0.005 ; t_tot=6*Nt_Sgro; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
signal =10000
stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
signal_trace=np.zeros(len(t))
signal_trace[stim_time_step:] = signal

# initializations
A_trace_orig_smalle=[A0]; R_trace_orig_smalle=[R0]
A_trace_orig_bige=[A0]; R_trace_orig_bige=[R0]

for i in range(len(t)-1):
    A_now=A_trace_orig_1[i]
    R_now=R_trace_orig_1[i]
    signal_now=signal_trace[i]
    
    A_next,R_next,r_now=Sgro_agent_smalle.update(dt,signal_now)
    A_trace_orig_smalle.append(A_next)
    R_trace_orig_smalle.append(R_next)
    
    A_next,R_next,r_now=Sgro_agent_bige.update(dt,signal_now)
    A_trace_orig_bige.append(A_next)
    R_trace_orig_bige.append(R_next)

A_trace_offset=1.5
A_trace_orig_smalle = np.array(A_trace_orig_smalle) # vectorize A_trace_orig
R_trace_orig_smalle = np.array(R_trace_orig_smalle)
A_trace_plot_smalle=(A_trace_orig_smalle+A_trace_offset)/Nh_Sgro

A_trace_orig_bige = np.array(A_trace_orig_bige) # vectorize A_trace_orig
R_trace_orig_bige = np.array(R_trace_orig_bige)
A_trace_plot_bige=(A_trace_orig_bige+A_trace_offset)/Nh_Sgro;

t_plot_Sgro = np.array(t)/(Nt_Sgro)
# nullclines
A_null = np.linspace(-2.5,2.5,num=200)
dAdt_null=A_null-1/3*A_null**3+a*np.log(1+signal/Kd)
dRdt_null=1/g*(A_null+c0)

#%% save output in npz files
np.savez('figS2_excitability_200422.npz',
         t_plot_Sgro = t_plot_Sgro,
         e_small = e_small, e_big = e_big,
         A_null = A_null, dRdt_null = dRdt_null,
         dAdt_null = dAdt_null, 
         
         A_trace_plot_smalle = A_trace_plot_smalle,
         A_trace_orig_smalle  =A_trace_orig_smalle, 
         R_trace_orig_smalle = R_trace_orig_smalle,
         A_trace_plot_bige = A_trace_plot_bige,
         A_trace_orig_bige  =A_trace_orig_bige, 
         R_trace_orig_bige = R_trace_orig_bige,
         A_mesh = A_mesh, R_mesh = R_mesh,
         dA = dA, dR = dR)

#%% load data
from iocustom import import_npz
import_npz('figS2_excitability_200422.npz',globals())

#%% Different excitability- trace, nullcline & streamplot
#label_font_size = 15
#trace_width = 2
tick_font_size = 20

fig5 = plt.figure(figsize=(12,9))
grid = plt.GridSpec(3, 2, wspace=0.3, hspace=0.8)

ax0 = fig5.add_subplot(grid[0, 0])
ax0.plot(t_plot_Sgro, A_trace_plot_smalle, 'green',linewidth=trace_width)
ax0.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
ax0.set_title('Excitability: '+str(e_small),fontsize=label_font_size)
ax0.set_xlim([0,6]); ax0.set_ylim([-0.4,1.2])
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax1= fig5.add_subplot(grid[1:, 0])
ax1.plot( A_null,dRdt_null,'darkgrey',linewidth=trace_width)
ax1.plot( A_null,dAdt_null,'deepskyblue',linewidth=trace_width)
ax1.plot( A_trace_orig_smalle, R_trace_orig_smalle,'darkgreen',linewidth=trace_width)
ax1.set_ylim([-0.7,3])
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)

A_arr = np.arange(-2.5,3,0.5)
R_arr = np.arange(-0.7,3,0.2)
makelong = 1
A_mesh, R_mesh = np.meshgrid(A_arr , R_arr )
dA = A_mesh-(np.power(A_mesh,3))/3-R_mesh+a*np.log(1+signal/Kd); dA=makelong*dA
dR = e_small*(A_mesh-g* R_mesh+c0); dR=makelong*dR
#ax1.quiver(A_mesh, R_mesh,dA,dR, linewidths=0.1, edgecolors='k')
ax1.streamplot(A_mesh,R_mesh,dA, dR,density=1,color='forestgreen')
ax1.plot(A_trace_orig_smalle[0],R_trace_orig_smalle[0],'*',markersize = 12,color='springgreen')


ax00 = fig5.add_subplot(grid[0, 1])
ax00.plot(t_plot_Sgro, A_trace_plot_bige, 'green',linewidth=trace_width)
ax00.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
ax00.set_title('Excitability: '+str(e_big),fontsize=label_font_size)
ax00.set_xlim([0,6]); ax00.set_ylim([-0.4,1.2])
ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax2= fig5.add_subplot(grid[1:, 1])
ax2.plot(A_null,dRdt_null, 'darkgrey',linewidth=trace_width)
ax2.plot(A_null, dAdt_null,'deepskyblue',linewidth=trace_width)
ax2.plot( A_trace_orig_bige, R_trace_orig_bige,'green',linewidth=trace_width)
ax2.set_ylim([-0.7,3])
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax1.set_aspect('equal')
A_arr = np.arange(-2.5,3,0.5)
R_arr = np.arange(-0.7,3,0.2)
makelong = 1
A_mesh, R_mesh = np.meshgrid(A_arr , R_arr )
dA = A_mesh-(np.power(A_mesh,3))/3-R_mesh+a*np.log(1+signal/Kd); dA=makelong*dA
dR = e_big*(A_mesh-g* R_mesh+c0); dR=makelong*dR
ax2.quiver(A_mesh, R_mesh,dA,dR, linewidths=0.1, edgecolors='k')
ax2.streamplot(A_mesh, R_mesh,dA,dR, density=1,color='forestgreen')
ax2.plot(A_trace_orig_bige[0],R_trace_orig_bige[0],'*',markersize = 12,color='springgreen')

fig5.text(0.04, 0.95,'A',ha='center', va='center', fontsize=abcd_font_size)
fig5.text(0.04, 0.8, r'$cAMP_{i}$'+'\n(Activator)', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.5, 0.04,r'$cAMP_{i}$'+'(Activator)',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.04, 0.35, 'Inhibitor', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

fig5.text(0.5, 0.95,'Sgro 2015',ha='center', va='center',color=mycolors[5],fontsize=title_font_size)
plt.show()

##%% vector field
#from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_agent_3var
#
#dt=0.001 ; t_tot=6*Nt_Goldbeter; t=list(np.arange(0,t_tot,dt))
#signal_small=0.05
#stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
#signal_trace_small=np.zeros(len(t))
#signal_trace_small[stim_time_step:] = signal_small
#
## initializations
#p0=0.897;b0=0.342
#p_trace_small=[p0]; b_trace_small=[b0]; g_trace_small=[g0]
#
## adaptive spike
#Goldbeter3_agent=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
#for i in range(len(t)-1):
#    p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,signal_trace_small[i])
#    p_trace_small.append(p_next)
#    b_trace_small.append(b_next)
#    g_trace_small.append(g_next)
#p_trace_small = np.array(p_trace_small); # p_trace_1 = p_trace_1/np.amax(p_trace_1)    
#b_trace_small = np.array(b_trace_small); # b_trace_1 = b_trace_1/np.amax(b_trace_1)    
#
#p_null = np.linspace(0.5, 1.25,num=50)
#signal = signal_small
#f1 = (k1+k2*signal)/(1+signal)
#f2 = (k1*L1+k2*L2*c*signal)/(1+c*signal)
#Ysq = (p_null*signal/(1+ signal))**2
#PI = a0*(lamda*theta + e*Ysq)/(1 + theta*a0 + (1 + a0)*e*Ysq)
#
#dpdt_null_small = f2/(f1+f2)
#dbdt_null_small = q*sig*PI/(ki+kt)
#trace_width = 3
#
#fig5 = plt.figure(figsize=(6,9))
#grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.8)
#
#ax0 = fig5.add_subplot(grid[0, 0])
#ax0.plot(t_plot_Goldbeter, b_trace_small, 'green',linewidth=trace_width)
#ax0.axvline(x=5, ls='--', linewidth=trace_width, color=mycolors[6])
#ax0.set_title(r'$cAMP_{e}$ input: '+str(signal_small),fontsize=label_font_size)
#ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
#ax0.set_xlim([0,30])
#
#ax1= fig5.add_subplot(grid[1:, 0])
#
##ax1.axvline(x=dpdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
#ax1.axvline(x=dpdt_null_small, ls='-', linewidth=trace_width, color='dimgrey')
#
##ax1.plot(p_null,dbdt_null_no_stim, 'lightblue',linewidth=trace_width)
#ax1.plot(p_null, dbdt_null_small,'deepskyblue',linewidth=trace_width)
#ax1.plot( p_trace_small, b_trace_small,'green',linewidth= 4)
#ax1.plot(p_trace_small[0],b_trace_small[0],'*',markersize = 12,color='springgreen')
#ax1.set_xlim([0.5,1.25]); ax1.set_ylim([-1,10])
#ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
#
## vector field for cAMPe = 0
#p_arr = np.arange(0,2,0.02)
#b_arr = np.arange(-5,15,0.2); makelong = 1
#p_mesh, b_mesh = np.meshgrid(p_arr , b_arr )
#signal = signal_small
#f1 = (k1+k2*signal)/(1+signal)
#f2 = (k1*L1+k2*L2*c*signal)/(1+c*signal)
#Ysq = (p_mesh*signal/(1+ signal))**2 
#PI = a0*(lamda*theta + e*Ysq)/(1 + theta*a0 + (1 + a0)*e*Ysq)
#dp = -f1*p_mesh+f2*(1-p_mesh);
#db = q*sig*PI-(kt+ki)*b_mesh
#ax1.streamplot(p_mesh,b_mesh,dp, db,color='forestgreen')
#
#fig5.text(0.04, 0.95,'B',ha='center', va='center', fontsize=abcd_font_size)
#
#fig5.text(0.03, 0.8,  r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
#fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)
#
#fig5.text(0.5, 0.04,'Inhibitor',ha='center', va='center',fontsize=label_font_size)
#fig5.text(0.03, 0.35, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
#
#fig5.text(0.5, 0.95,'Martiel 1987',ha='center', va='center',color=mycolors[0],fontsize=title_font_size)
#plt.show()
#%% Goldbeter 1987 nullclines
from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_agent_3var

# Table 2 parameters
k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10; L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01; theta=0.01
e=  1 
q=4000
sig= 0.6
v=12; k= 4 # k prime in the paper
ki=1.7
kt=0.9
kc=5.4
h=5

Goldbeter3AgentParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}
p0=0.9; a0=3; b0=0.2; g0=0

dt=0.001; t_tot=6*Nt_Goldbeter; t=list(np.arange(0,t_tot,dt))
# define extracellular stim trace
#constant_signal=1 
signal_1=1
signal_10k = 10000
stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
signal_trace_1=np.zeros(len(t))
signal_trace_1[stim_time_step:] = signal_1
signal_trace_10k=np.zeros(len(t))
signal_trace_10k[stim_time_step:] = signal_10k

# initializations
p_trace_1=[p0]; b_trace_1=[b0]; g_trace_1=[g0]
p_trace_10k=[p0]; b_trace_10k=[b0]; g_trace_10k=[g0]

# adaptive spike
# Fig 2, 4 variable model, autonomous oscillation
Goldbeter3_agent=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
for i in range(len(t)-1):
    p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,signal_trace_1[i])
    p_trace_1.append(p_next)
    b_trace_1.append(b_next)
    g_trace_1.append(g_next)
p_trace_1 = np.array(p_trace_1);    
b_trace_1 = np.array(b_trace_1); 
b_trace_1_hnorm = b_trace_1/Nh_Goldbeter

# sustained oscillations
Goldbeter3_agent=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
for i in range(len(t)-1):
    p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,signal_trace_10k[i])
    p_trace_10k.append(p_next)
    b_trace_10k.append(b_next)
    g_trace_10k.append(g_next)
p_trace_10k = np.array(p_trace_10k);    
b_trace_10k = np.array(b_trace_10k); 
b_trace_10k_hnorm = b_trace_10k/Nh_Goldbeter

t_plot_Goldbeter = np.array(t)/(Nt_Goldbeter)

# nullclines
p_null = np.linspace(-2, 6,num=200)
signal = 0
f1 = (k1+k2*signal)/(1+signal)
f2 = (k1*L1+k2*L2*c*signal)/(1+c*signal)
Ysq = (p_null*signal/(1+ signal))**2 
PI = a0*(lamda*theta + e*Ysq)/(1 + theta*a0 + (1 + a0)*e*Ysq)

dpdt_null_no_stim = f2/(f1+f2)
dbdt_null_no_stim = q*sig*PI/(ki+kt)

signal = signal_1
f1 = (k1+k2*signal)/(1+signal)
f2 = (k1*L1+k2*L2*c*signal)/(1+c*signal)
Ysq = (p_null*signal/(1+ signal))**2
PI = a0*(lamda*theta + e*Ysq)/(1 + theta*a0 + (1 + a0)*e*Ysq)

dpdt_null_1 = f2/(f1+f2)
dbdt_null_1 = q*sig*PI/(ki+kt)

signal = signal_10k
f1 = (k1+k2*signal)/(1+signal)
f2 = (k1*L1+k2*L2*c*signal)/(1+c*signal)
Ysq = (p_null*signal/(1+ signal))**2
PI = a0*(lamda*theta + e*Ysq)/(1 + theta*a0 + (1 + a0)*e*Ysq)

dpdt_null_10k = f2/(f1+f2)
dbdt_null_10k = q*sig*PI/(ki+kt)

# adaptive spike & oscillation- trace and nullcline

#label_font_size = 15
#trace_width = 3


#%% 
fig5 = plt.figure(figsize=(12,9))
grid = plt.GridSpec(3, 2, wspace=0.3, hspace=0.8)

ax0 = fig5.add_subplot(grid[0, 0])
ax0.plot(t_plot_Goldbeter, b_trace_1_hnorm, 'green',linewidth=trace_width)
ax0.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
ax0.set_title(r'$cAMP_{e}$ input: '+str(signal_1),fontsize=label_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_xlim([0,6]); ax0.set_ylim([-0.2,2])
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax1= fig5.add_subplot(grid[1:, 0])

ax1.axvline(x=dpdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax1.axvline(x=dpdt_null_1, ls='-', linewidth=trace_width, color='dimgrey')

ax1.plot(p_null,dbdt_null_no_stim, 'lightblue',linewidth=trace_width)
ax1.plot(p_null, dbdt_null_1,'deepskyblue',linewidth=trace_width)
ax1.plot( p_trace_1, b_trace_1,'green',linewidth= 4)
ax1.plot(p_trace_1[0],b_trace_1[0],'*',markersize = 12,color='springgreen')
ax1.set_xlim([-0.5,1.5]); ax1.set_ylim([-20,700])
# ax1.set_ylim([0,0.5]); 
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax00 = fig5.add_subplot(grid[0, 1])
ax00.plot(t_plot_Goldbeter, b_trace_10k_hnorm, 'green',linewidth=trace_width)
ax00.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
ax00.set_title(r'$cAMP_{e}$ input: '+str(signal_10k),fontsize=label_font_size)
ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax00.set_xlim([0,6]); ax00.set_ylim([-0.2,2])
ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax2= fig5.add_subplot(grid[1:, 1])

ax2.axvline(x=dpdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax2.axvline(x=dpdt_null_10k, ls='-', linewidth=trace_width, color='dimgrey')

ax2.plot(p_null,dbdt_null_no_stim, 'lightblue',linewidth=trace_width)
ax2.plot(p_null, dbdt_null_10k,'deepskyblue',linewidth=trace_width)
ax2.plot( p_trace_10k, b_trace_10k,'green',linewidth= 4)
ax2.plot(p_trace_10k[0],b_trace_10k[0],'*',markersize = 12,color='springgreen')
ax2.set_xlim([-0.5,1.5]); ax2.set_ylim([-20,700])
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)

fig5.text(0.04, 0.95,'A',ha='center', va='center', fontsize=abcd_font_size)

fig5.text(0.04, 0.8,  r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)

fig5.text(0.5, 0.04,'Inhibitor',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.04, 0.35, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

fig5.text(0.5, 0.95,'Martiel 1987',ha='center', va='center',color=mycolors[0],fontsize=title_font_size)
plt.show()

#%% Make videos
fig5 = plt.figure(figsize=(12,9))
grid = plt.GridSpec(3, 2, wspace=0.3, hspace=0.8)
fig5.text(0.04, 0.95,'A',ha='center', va='center', fontsize=abcd_font_size)

fig5.text(0.04, 0.95,'A',ha='center', va='center', fontsize=abcd_font_size)

fig5.text(0.04, 0.8,  r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)

fig5.text(0.5, 0.04,'Inhibitor',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.04, 0.35, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

fig5.text(0.5, 0.95,'Martiel 1987',ha='center', va='center',color=mycolors[0],fontsize=title_font_size)
ax0 = fig5.add_subplot(grid[0, 0])
ax1= fig5.add_subplot(grid[1:, 0])
ax00 = fig5.add_subplot(grid[0, 1])
ax2= fig5.add_subplot(grid[1:, 1])

camera = Camera(fig5)
down_sample_frame =  350

for huhu in range(int(len(t_plot_Sgro)/down_sample_frame)): 
    i = huhu*down_sample_frame
    ax0.plot(t_plot_Goldbeter[0:i], b_trace_1_hnorm[0:i], 'green',linewidth=trace_width)
    ax0.plot(t_plot_Goldbeter[i], b_trace_1_hnorm[i],'o',markersize = 8,color='springgreen')
    ax0.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
    ax0.set_title(r'$cAMP_{e}$ input: '+str(signal_1),fontsize=label_font_size)
    ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax0.set_xlim([0,6]); ax0.set_ylim([-0.2,2])
    ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    if i < stim_time_step:
        ax1.axvline(x=dpdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
        ax1.plot(p_null,dbdt_null_no_stim, 'lightblue',linewidth=trace_width)
    else:
        ax1.axvline(x=dpdt_null_1, ls='-', linewidth=trace_width, color='dimgrey')
        ax1.plot(p_null, dbdt_null_1,'deepskyblue',linewidth=trace_width)
    ax1.plot( p_trace_1[:i], b_trace_1[:i],'green',linewidth= 4)
    ax1.plot(p_trace_1[i],b_trace_1[i],'o',markersize = 8,color='springgreen')
    ax1.set_xlim([-0.5,1.5]); ax1.set_ylim([-20,700])
    # ax1.set_ylim([0,0.5]); 
    ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    ax00.plot(t_plot_Goldbeter[0:i], b_trace_10k_hnorm[0:i], 'green',linewidth=trace_width)
    ax00.plot(t_plot_Goldbeter[i], b_trace_10k_hnorm[i],'o',markersize = 8,color='springgreen')
    ax00.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
    ax00.set_title(r'$cAMP_{e}$ input: '+str(signal_10k),fontsize=label_font_size)
    ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax00.set_xlim([0,6]); ax00.set_ylim([-0.2,2])
    ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    if i < stim_time_step:
        ax2.axvline(x=dpdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
        ax2.plot(p_null,dbdt_null_no_stim, 'lightblue',linewidth=trace_width)
    else:
        ax2.axvline(x=dpdt_null_10k, ls='-', linewidth=trace_width, color='dimgrey')
        ax2.plot(p_null, dbdt_null_10k,'deepskyblue',linewidth=trace_width)
    ax2.plot( p_trace_10k[0:i], b_trace_10k[0:i],'green',linewidth= 4)
    ax2.plot(p_trace_10k[i],b_trace_10k[i],'o',markersize = 8,color='springgreen')
    ax2.set_xlim([-0.5,1.5]); ax2.set_ylim([-20,700])
    ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    camera.snap()     
ani = camera.animate()  
#%%  
OutPath = 'C:\\Users\\ellin\\Dropbox\\Prospectus\\Presentation\\nullcline_videos\\'
ani_name = 'Martiel_bif'
ani.save(OutPath + ani_name + '.mp4') 
#%% Kamino nullclines, adaptive spike and oscillations

from Kamino2017_agent_and_pop_FUN import Kamino2017_agent

tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
AgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}

dt=0.001 ; t_tot=6*Nt_Kamino; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
#constant_signal=1 
signal_1=1
signal_10k = 10000
stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
signal_trace_1=np.zeros(len(t))
signal_trace_1[stim_time_step:] = signal_1
signal_trace_10k=np.zeros(len(t))
signal_trace_10k[stim_time_step:] = signal_10k

# initializations
x0=0.01; y0=0.05; z0=0.005
x_trace_1=[x0]; y_trace_1=[y0]
x_trace_10k=[x0]; y_trace_10k=[y0]

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

# sustained oscillations
Kamino_agent=Kamino2017_agent([x0,y0,z0],AgentParam)
for i in range(len(t)-1):
    x_now=x_trace_10k[i]
    y_now=y_trace_10k[i]
    x_next,y_next,z_next= Kamino_agent.update(dt, signal_trace_10k[i])
    x_trace_10k.append(x_next)
    y_trace_10k.append(y_next)

x_trace_10k = np.array(x_trace_10k)
y_trace_10k = np.array(y_trace_10k) 
y_trace_10k_hnorm = (y_trace_10k-Nh_Kamino_offset)/Nh_Kamino

t_plot_Kamino = np.array(t)/(Nt_Kamino)

# nullclines
x_null_short = np.linspace(-0.5,2,num=100)
x_null_long = np.linspace(-1000,11000,num=100)
dxdt_null_no_stim = 0 + delta
dxdt_null_1 = signal_1 + delta
dxdt_null_10k = signal_10k + delta 
dydt_null_no_stim_short = (0+delta)**n/((0+delta)**n+(np.power(K*x_null_short,n)))
dydt_null_no_stim_long = (0+delta)**n/((0+delta)**n+(np.power(K*x_null_long,n)))
dydt_null_1 = (signal_1+delta)**n/((signal_1+delta)**n+(np.power(K*x_null_short,n)))
dydt_null_10k = (signal_10k+delta)**n/((signal_10k+delta)**n+(np.power(K*x_null_long,n)))


#%%adaptive spike & oscillation- trace and nullcline
abcd_font_size = 28
label_font_size=24
title_font_size = 26
sublabel_font_size = 22
trace_width=3
tick_font_size=20

fig5 = plt.figure(figsize=(12,9))
grid = plt.GridSpec(3, 2, wspace=0.3, hspace=0.8)

ax0 = fig5.add_subplot(grid[0, 0])
ax0.plot(t_plot_Kamino, y_trace_1_hnorm, 'green',linewidth=trace_width)
ax0.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
ax0.set_title(r'$cAMP_{e}$ input: '+str(signal_1),fontsize=label_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_xlim([0,6]); ax0.set_ylim([-0.2,1.2])
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax1= fig5.add_subplot(grid[1:, 0])
ax1.axvline(x=dxdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax1.axvline(x=dxdt_null_1, ls='-', linewidth=trace_width, color='dimgrey')
ax1.plot(x_null_short,dydt_null_no_stim_short,'lightblue',linewidth=trace_width)
ax1.plot(x_null_short, dydt_null_1,'deepskyblue',linewidth=trace_width)

ax1.plot( x_trace_1, y_trace_1,'green',linewidth=trace_width)
ax1.plot( x_trace_1[0], y_trace_1[0],'*',markersize = 12,color='springgreen')
ax1.set_xlim([-0.5,2])
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
#ax1.set_aspect('equal')


ax00 = fig5.add_subplot(grid[0, 1])
ax00.plot(t_plot_Kamino, y_trace_10k_hnorm, 'green',linewidth=trace_width)
ax00.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
ax00.set_title(r'$cAMP_{e}$ input: '+str(signal_10k),fontsize=label_font_size)
ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax00.set_xlim([0,6]); ax00.set_ylim([-0.2,1.2])
ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax2= fig5.add_subplot(grid[1:, 1])
ax2.axvline(x=dxdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax2.axvline(x=dxdt_null_10k, ls='-', linewidth=trace_width, color='dimgrey')
ax2.plot(x_null_long,dydt_null_no_stim_long,'lightblue',linewidth=trace_width)
ax2.plot(x_null_long, dydt_null_10k,'deepskyblue',linewidth=trace_width)
ax2.plot( x_trace_10k, y_trace_10k,'green',linewidth=trace_width)
ax2.plot( x_trace_10k[0], y_trace_10k[0],'*',markersize = 12,color='springgreen')
ax2.set_xlim([-1000,11000])
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
#ax2.set_xlabel( 'Activator' ,fontsize=label_font_size)
#ax2.set_ylabel( 'Inhibitor' ,fontsize=label_font_size)
#ax2.set_aspect('equal')

#cAMP_ext = 10000
#A_arr = np.arange(-2.5,2.5,0.2)
#R_arr = np.arange(-1,2.5,0.2)
#A_mesh, R_mesh = np.meshgrid(A_arr , R_arr )
#dA = A_mesh-(np.power(A_mesh,3))/3-R_mesh+a*np.log(1+cAMP_ext/Kd)
#dR = e*(A_mesh-g* R_mesh+c0)
#ax2.quiver(A_mesh, R_mesh,dA,dR, linewidths=0.1, edgecolors='k')
fig5.text(0.04, 0.95,'C',ha='center', va='center', fontsize=abcd_font_size)

fig5.text(0.04, 0.8,  r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.5, 0.04,'Inhibitor',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.04, 0.35, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

fig5.text(0.5, 0.95,'Kamino 2017',ha='center', va='center',color=mycolors[7],fontsize=title_font_size)
plt.show()

#%% nullclines make video
fig5 = plt.figure(figsize=(12,9))
grid = plt.GridSpec(3, 2, wspace=0.3, hspace=0.8)
fig5.text(0.04, 0.95,'C',ha='center', va='center', fontsize=abcd_font_size)

fig5.text(0.04, 0.8,  r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.5, 0.04,'Inhibitor',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.04, 0.35, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.95,'Kamino 2017',ha='center', va='center',color=mycolors[7],fontsize=title_font_size)
ax0 = fig5.add_subplot(grid[0, 0])
ax1= fig5.add_subplot(grid[1:, 0])
ax00 = fig5.add_subplot(grid[0, 1])
ax2= fig5.add_subplot(grid[1:, 1])

camera = Camera(fig5)
down_sample_frame =  350

for huhu in range(int(len(t_plot_Kamino)/down_sample_frame)): 
    i = huhu*down_sample_frame
    ax0.plot(t_plot_Kamino[0:i], y_trace_1_hnorm[0:i], 'green',linewidth=trace_width)
    ax0.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
    ax0.set_title(r'$cAMP_{e}$ input: '+str(signal_1),fontsize=label_font_size)
    ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax0.set_xlim([0,6]); ax0.set_ylim([-0.2,1.2])
    ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    ax1= fig5.add_subplot(grid[1:, 0])
    if i < stim_time_step:
        ax1.axvline(x=dxdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
        ax1.plot(x_null_short,dydt_null_no_stim_short,'lightblue',linewidth=trace_width)
    else:
        ax1.axvline(x=dxdt_null_1, ls='-', linewidth=trace_width, color='dimgrey')
        ax1.plot(x_null_short, dydt_null_1,'deepskyblue',linewidth=trace_width)
    
    ax1.plot( x_trace_1[0:i], y_trace_1[0:i],'green',linewidth=trace_width)
    ax1.plot( x_trace_1[i], y_trace_1[i],'o',markersize = 8,color='springgreen')
    ax1.set_xlim([-0.5,2])
    ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
    #ax1.set_aspect('equal')
    
    
    ax00 = fig5.add_subplot(grid[0, 1])
    ax00.plot(t_plot_Kamino[0:i], y_trace_10k_hnorm[0:i], 'green',linewidth=trace_width)
    ax00.axvline(x=1, ls='--', linewidth=trace_width, color=mycolors[6])
    ax00.set_title(r'$cAMP_{e}$ input: '+str(signal_10k),fontsize=label_font_size)
    ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax00.set_xlim([0,6]); ax00.set_ylim([-0.2,1.2])
    ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    ax2= fig5.add_subplot(grid[1:, 1])
    if i < stim_time_step:
        ax2.axvline(x=dxdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
        ax2.plot(x_null_long,dydt_null_no_stim_long,'lightblue',linewidth=trace_width)
    else:
        ax2.axvline(x=dxdt_null_10k, ls='-', linewidth=trace_width, color='dimgrey')
        ax2.plot(x_null_long, dydt_null_10k,'deepskyblue',linewidth=trace_width)
    ax2.plot( x_trace_10k[0:i], y_trace_10k[0:i],'green',linewidth=trace_width)
    ax2.plot( x_trace_10k[i], y_trace_10k[i],'o',markersize = 8,color='springgreen')
    ax2.set_xlim([-1000,11000])
    ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
    camera.snap()       
ani = camera.animate()
 #%%
OutPath = 'C:\\Users\\ellin\\Dropbox\\Prospectus\\Presentation\\nullcline_videos\\'
ani_name = 'Kamino_bif'
ani.save(OutPath + ani_name + '.mp4') 
     


#%% save figure S1 output in npz files
np.savez('figS1_nullcline_200422.npz',
         t_plot_Sgro = t_plot_Sgro,
         signal_1 = signal_1, signal_10k = signal_10k,
         A_trace_plot_1 = A_trace_plot_1,
         A_trace_plot_10k = A_trace_plot_10k,
         A_null = A_null, dRdt_null = dRdt_null,
         dAdt_null_no_stim = dAdt_null_no_stim, 
         dAdt_null_1 = dAdt_null_1, dAdt_null_10k = dAdt_null_10k,
         A_trace_orig_1 = A_trace_orig_1 , R_trace_orig_1 = R_trace_orig_1,
         A_trace_orig_10k = A_trace_orig_10k, R_trace_orig_10k = R_trace_orig_10k,
         
         t_plot_Goldbeter = t_plot_Goldbeter, 
         b_trace_1_hnorm = b_trace_1_hnorm, b_trace_10k_hnorm = b_trace_10k_hnorm,
         dpdt_null_no_stim = dpdt_null_no_stim, 
         dpdt_null_1 = dpdt_null_1, dpdt_null_10k = dpdt_null_10k,
         p_null = p_null, dbdt_null_no_stim = dbdt_null_no_stim,
         dbdt_null_1 = dbdt_null_1, dbdt_null_10k = dbdt_null_10k, 
         p_trace_1 = p_trace_1, b_trace_1 =  b_trace_1,
         p_trace_10k = p_trace_10k, b_trace_10k = b_trace_10k,
         
         t_plot_Kamino = t_plot_Kamino, 
         y_trace_1_hnorm = y_trace_1_hnorm, y_trace_10k_hnorm = y_trace_10k_hnorm,
         dxdt_null_no_stim = dxdt_null_no_stim,
         dxdt_null_1 = dxdt_null_1, dxdt_null_10k = dxdt_null_10k,
         x_null_short = x_null_short, dydt_null_no_stim_short = dydt_null_no_stim_short,
         dydt_null_1 = dydt_null_1, dydt_null_10k = dydt_null_10k,
         x_trace_1 = x_trace_1, y_trace_1= y_trace_1,
         x_null_long = x_null_long, dydt_null_no_stim_long = dydt_null_no_stim_long,
         x_trace_10k = x_trace_10k, y_trace_10k = y_trace_10k)


#%% Goldbeter, vector field, with ramp cAMP input
from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_agent_3var
# Table 2 parameters
k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10; L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01; theta=0.01
e=  1 
q=4000
sig= 0.6
v=12; k= 4 # k prime in the paper
ki=1.7
kt=0.9
kc=5.4
h=5
Goldbeter3AgentParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}
# initializations
p0=0.897;b0=0.342; g0 = 0; a0=3
p_trace=[p0]; b_trace=[b0]; g_trace=[g0]

Goldbeter3_agent=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
dt=0.001
# define extracellular stim trace
signal = 1
Exp_time = Sgro2015Figure3excel["Ramp Input (min Time)"]# np.linspace(0,20,num=1000) # t step 0.02
Exp_time = Exp_time[~np.isnan(Exp_time)]/4*Nt_Goldbeter
RampInput_Exp = Sgro2015Figure3excel["Ramp Input (nM cAMP)"]
RampInput_Exp = RampInput_Exp[~np.isnan(RampInput_Exp)]
RampInput_Exp = RampInput_Exp[math.floor(len(RampInput_Exp)*4/8): math.floor(len(RampInput_Exp)*7/8)]
Exp_time = Exp_time[0:len(RampInput_Exp)]

Goldbeter_time = np.arange(0,7.5*Nt_Goldbeter,dt)
RampInput_Goldbeter= np.interp(Goldbeter_time, Exp_time,RampInput_Exp) * signal

for i in range(len(Goldbeter_time)-1):
    signal_now = RampInput_Goldbeter[i]
    p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,signal_now)
    p_trace.append(p_next)
    b_trace.append(b_next)
    g_trace.append(g_next)
p_trace = np.array(p_trace);  
b_trace = np.array(b_trace);  
b_trace_hnorm = b_trace/Nh_Goldbeter
t_plot_Goldbeter = np.array(Goldbeter_time)/(Nt_Goldbeter)

# nullclines
p_null = np.linspace(-2, 6,num=200)

signal = 0
f1 = (k1+k2*signal)/(1+signal)
f2 = (k1*L1+k2*L2*c*signal)/(1+c*signal)
Ysq = (p_null*signal/(1+ signal))**2 
PI = a0*(lamda*theta + e*Ysq)/(1 + theta*a0 + (1 + a0)*e*Ysq)

dpdt_null_no_stim = f2/(f1+f2)
dbdt_null_no_stim = q*sig*PI/(ki+kt)

signal = 1
f1 = (k1+k2*signal)/(1+signal)
f2 = (k1*L1+k2*L2*c*signal)/(1+c*signal)
Ysq = (p_null*signal/(1+ signal))**2
PI = a0*(lamda*theta + e*Ysq)/(1 + theta*a0 + (1 + a0)*e*Ysq)

dpdt_null = f2/(f1+f2)
dbdt_null = q*sig*PI/(ki+kt)

trace_width = 3
         
         
#%%
# fig5 = plt.figure(figsize=(8,9))
# grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.8)
abcd_font_size = 28
label_font_size=20
title_font_size = 22
sublabel_font_size = 20
trace_width=3
tick_font_size=18
fig5 = plt.figure(figsize=(8,7))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.55)

ax0 = fig5.add_subplot(grid[0, 0],xticks=[0,2.5,5,7.5])
ax0.plot(t_plot_Goldbeter, b_trace_hnorm, 'green',linewidth=trace_width)
ax0.axvspan(2.5, 7.5, alpha=0.2, color='b');
ax0.set_title(r'$cAMP_{e}$ ramp input',fontsize=label_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_xlim([0,7.5]);  ax0.set_ylim([-0.03,0.18])
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax1= fig5.add_subplot(grid[1:, 0])

#ax1.axvline(x=dpdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax1.axvline(x=dpdt_null, ls='-', linewidth=trace_width, color='dimgrey')
#ax1.plot(p_null,dbdt_null_no_stim, 'lightblue',linewidth=trace_width)
ax1.plot(p_null, dbdt_null,'deepskyblue',linewidth=trace_width)

ax1.axvline(x=dpdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax1.axvline(x=dpdt_null, ls='-', linewidth=trace_width, color='dimgrey')

ax1.plot(p_null,dbdt_null_no_stim, 'lightblue',linewidth=trace_width)
ax1.plot(p_null, dbdt_null,'deepskyblue',linewidth=trace_width)

ax1.plot( p_trace, b_trace,'green',linewidth= 4)
ax1.plot(p_trace[0],b_trace[0],'*',markersize = 12,color='springgreen')
ax1.set_xlim([-0.5,1.5]); ax1.set_ylim([-30,600])
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)

# vector field for cAMPe = 1
p_arr = np.arange(-1,2,0.02)
b_arr = np.arange(-50,700,10); makelong = 1
p_mesh, b_mesh = np.meshgrid(p_arr , b_arr )
signal = 1
f1 = (k1+k2*signal)/(1+signal)
f2 = (k1*L1+k2*L2*c*signal)/(1+c*signal)
Ysq = (p_mesh*signal/(1+ signal))**2 
PI = a0*(lamda*theta + e*Ysq)/(1 + theta*a0 + (1 + a0)*e*Ysq)
dp = -f1*p_mesh+f2*(1-p_mesh);
db = q*sig*PI-(kt+ki)*b_mesh
ax1.streamplot(p_mesh,b_mesh,dp, db,color='forestgreen')

fig5.text(0.04, 0.95,'A',ha='center', va='center', fontsize=abcd_font_size)

fig5.text(0.03, 0.8,  r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)

fig5.text(0.5, 0.04,'Inhibitor',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.03, 0.35, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

fig5.text(0.5, 0.95,'Martiel 1987',ha='center', va='center',color=mycolors[0],fontsize=title_font_size)
plt.show()

#%% make movie
fig5 = plt.figure(figsize=(8,7))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.55)

fig5.text(0.04, 0.95,'A',ha='center', va='center', fontsize=abcd_font_size)
fig5.text(0.03, 0.8,  r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.5, 0.04,'Inhibitor',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.03, 0.35, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.95,'Martiel 1987',ha='center', va='center',color=mycolors[0],fontsize=title_font_size)

ax0 = fig5.add_subplot(grid[0, 0],xticks=[0,2.5,5,7.5])
ax1= fig5.add_subplot(grid[1:, 0])

camera = Camera(fig5)
down_sample_frame =  500

for huhu in range(int(len(t_plot_Goldbeter)/down_sample_frame)): 
    i = huhu*down_sample_frame
    
    ax0.plot(t_plot_Goldbeter[0:i], b_trace_hnorm[0:i], 'green',linewidth=trace_width)
    ax0.axvspan(2.5, 7.5, alpha=0.2, color='b');
    ax0.set_title(r'$cAMP_{e}$ ramp input',fontsize=label_font_size)
    ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax0.set_xlim([0,7.5]);  ax0.set_ylim([-0.03,0.18])
    ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    signal = RampInput_Goldbeter[i]
    f1 = (k1+k2*signal)/(1+signal)
    f2 = (k1*L1+k2*L2*c*signal)/(1+c*signal)
    Ysq = (p_null*signal/(1+ signal))**2
    PI = a0*(lamda*theta + e*Ysq)/(1 + theta*a0 + (1 + a0)*e*Ysq)
    dpdt_null = f2/(f1+f2)
    dbdt_null = q*sig*PI/(ki+kt)
    
    ax1.axvline(x=dpdt_null, ls='-', linewidth=trace_width, color='dimgrey')
    ax1.plot(p_null, dbdt_null,'deepskyblue',linewidth=trace_width)
    
    ax1.plot( p_trace[0:i], b_trace[0:i],'green',linewidth= 4)
    ax1.plot(p_trace[i],b_trace[i],'o',markersize = 8,color='springgreen')
    ax1.set_xlim([-0.5,1.5]); ax1.set_ylim([-50,600])
    ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    # # vector field 
    # p_arr = np.arange(-1,2,0.02)
    # b_arr = np.arange(-50,700,10); makelong = 1
    # p_mesh, b_mesh = np.meshgrid(p_arr , b_arr )

    # f1 = (k1+k2*signal)/(1+signal)
    # f2 = (k1*L1+k2*L2*c*signal)/(1+c*signal)
    # Ysq = (p_mesh*signal/(1+ signal))**2 
    # PI = a0*(lamda*theta + e*Ysq)/(1 + theta*a0 + (1 + a0)*e*Ysq)
    # dp = -f1*p_mesh+f2*(1-p_mesh);
    # db = q*sig*PI-(kt+ki)*b_mesh
    # ax1.streamplot(p_mesh,b_mesh,dp, db,color='forestgreen')

    camera.snap()    
ani = camera.animate()
   #%%
OutPath = 'C:\\Users\\ellin\\Dropbox\\Prospectus\\Presentation\\nullcline_videos\\'
ani_name = 'Goldbeter_rap_wo_streams'
ani.save(OutPath + ani_name + '.mp4') 
#%% Sgro 2015, vector field with ramp cAMP input
from Sgro2015_agent_and_pop_FUN import Sgro2015_agent

e=0.1; 
tauA=0.09; tauR=tauA/e; g=0.5; c0 = 1.2
a= 0.058; Kd = 1e-5
Na=3.5;  # normalization factor of A
A0=-1.5; R0=-0.5
Nt_Sgro = 27
dt=0.001 
# Create agent with small or big e
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)

# define extracellular stim trace
signal = 1
Exp_time = Sgro2015Figure3excel["Ramp Input (min Time)"]# np.linspace(0,20,num=1000) # t step 0.02
Exp_time = Exp_time[~np.isnan(Exp_time)]/4*Nt_Sgro
RampInput_Exp = Sgro2015Figure3excel["Ramp Input (nM cAMP)"]
RampInput_Exp = RampInput_Exp[~np.isnan(RampInput_Exp)]
RampInput_Exp = RampInput_Exp[math.floor(len(RampInput_Exp)*4/8): math.floor(len(RampInput_Exp)*7/8)]
Exp_time = Exp_time[0:len(RampInput_Exp)]

Sgro_time = np.arange(0,7.5*Nt_Sgro,dt)
RampInput_Sgro= np.interp(Sgro_time, Exp_time,RampInput_Exp) * signal

#plt.figure()
#plt.plot(Sgro_time/27,RampInput_Sgro)
#plt.show()

# initializations
A_trace_orig=[A0]; R_trace_orig=[R0]

for i in range(len(Sgro_time)-1):
    A_now=A_trace_orig[i]
    R_now=R_trace_orig[i]
    signal_now=RampInput_Sgro[i]
    
    A_next,R_next,r_now=Sgro_agent.update(dt,signal_now)
    A_trace_orig.append(A_next)
    R_trace_orig.append(R_next)
    

A_trace_offset=1.5
A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
R_trace_orig = np.array(R_trace_orig)
A_trace_plot=(A_trace_orig+A_trace_offset)/Na
t_plot_Sgro = np.array(Sgro_time)/(Nt_Sgro)
# nullclines
A_null = np.linspace(-2.5,2.5,num=200)
dAdt_null=A_null-1/3*A_null**3+a*np.log(1+signal/Kd)
dAdt_null_no_stim=A_null-1/3*(A_null**3)
dRdt_null=1/g*(A_null+c0)
#%%
# np.savetxt('Ramp_Sgro_nullcline_OUT_1204.txt', (t_plot_Sgro, A_trace_plot,A_trace_orig,R_trace_orig))
# read out the saved traces
[t_plot_Sgro, A_trace_plot,A_trace_orig,R_trace_orig] = np.loadtxt('Ramp_Sgro_nullcline_OUT_1204.txt', dtype=float)

#%% plot out outputs
trace_width = 3
fig5 = plt.figure(figsize=(8,9))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.8)

ax0 = fig5.add_subplot(grid[0, 0],xticks=[0,2.5,5,7.5])
ax0.plot(t_plot_Sgro, A_trace_plot, 'green',linewidth=trace_width)
# ax0.axvline(x=5, ls='--', linewidth=trace_width, color='dimgrey')
ax0.axvspan(2.5, 7.5, alpha=0.2, color='b');
ax0.set_title(r'$cAMP_{e}$ ramp input',fontsize=label_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_xlim([0,7.5]); ax0.set_ylim([-0.2,1.2])
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax1= fig5.add_subplot(grid[1:, 0])

#ax1.plot(p_null,dbdt_null_no_stim, 'lightblue',linewidth=trace_width)
ax1.plot(A_null, dRdt_null,'dimgrey',linewidth=trace_width)
ax1.plot(A_null, dAdt_null,'deepskyblue',linewidth=trace_width)
ax1.plot(A_null,dAdt_null_no_stim,'lightblue',linewidth=trace_width)
ax1.plot( A_trace_orig, R_trace_orig,'green',linewidth= 4)
ax1.plot(A_trace_orig[0],R_trace_orig[0],'*',markersize = 12,color='springgreen')
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_xlim([-2.5,2.5]); ax1.set_ylim([-1,3])
# vector field for tiny cAMPe inputs 

A_arr = np.arange(-2.5,3,0.5)
R_arr = np.arange(-1.2,3,0.5)
makelong = 1
A_mesh, R_mesh = np.meshgrid(A_arr , R_arr )
dA = A_mesh-(np.power(A_mesh,3))/3-R_mesh+a*np.log(1+signal/Kd); dA=makelong*dA
dR = e*(A_mesh-g* R_mesh+c0); dR=makelong*dR

ax1.streamplot(A_mesh,R_mesh,dA, dR,color='forestgreen')
fig5.text(0.04, 0.95,'B',ha='center', va='center', fontsize=abcd_font_size)

fig5.text(0.03, 0.8,  r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)

fig5.text(0.5, 0.04,'Inhibitor',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.03, 0.35, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

fig5.text(0.5, 0.95,'Sgro 2015',ha='center', va='center',color=mycolors[5],fontsize=title_font_size)
plt.show()

#%% Make video Sgro
fig5 = plt.figure(figsize=(8,7))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.55)

fig5.text(0.04, 0.95,'B',ha='center', va='center', fontsize=abcd_font_size)
fig5.text(0.03, 0.8,  r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.5, 0.04,r'$cAMP_{i}$',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.03, 0.35, 'Inhibitor', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.95,'Sgro 2015',ha='center', va='center',color=mycolors[5],fontsize=title_font_size)

ax0 = fig5.add_subplot(grid[0, 0],xticks=[0,2.5,5,7.5])
ax1= fig5.add_subplot(grid[1:, 0])

camera = Camera(fig5)
down_sample_frame =  500

for huhu in range(int(len(t_plot_Sgro)/down_sample_frame)): 
    i = huhu*down_sample_frame
    ax0.plot(t_plot_Sgro[0:i], A_trace_plot[0:i], 'green',linewidth=trace_width)
    ax0.axvspan(2.5, 7.5, alpha=0.2, color='b');
    ax0.set_title(r'$cAMP_{e}$ ramp input',fontsize=label_font_size)
    ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax0.set_xlim([0,7.5]); ax0.set_ylim([-0.2,1.2])
    ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    signal = RampInput_Sgro[i]
    A_null = np.linspace(-2.5,2.5,num=200)
    dAdt_null=A_null-1/3*A_null**3+a*np.log(1+signal/Kd)
    dRdt_null=1/g*(A_null+c0)
    
    ax1.plot(A_null, dRdt_null,'dimgrey',linewidth=trace_width)
    ax1.plot(A_null, dAdt_null,'deepskyblue',linewidth=trace_width)
    
    ax1.plot( A_trace_orig[0:i], R_trace_orig[0:i],'green',linewidth= 4)
    ax1.plot(A_trace_orig[i],R_trace_orig[i],'o',markersize = 8,color='springgreen')
    ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax1.set_xlim([-2.5,2.5]); ax1.set_ylim([-1,3])
    
    # vector field for tiny cAMPe inputs 
    A_arr = np.arange(-2.5,3,0.5)
    R_arr = np.arange(-1.2,3,0.5)
    makelong = 1
    A_mesh, R_mesh = np.meshgrid(A_arr , R_arr )
    dA = A_mesh-(np.power(A_mesh,3))/3-R_mesh+a*np.log(1+signal/Kd); dA=makelong*dA
    dR = e*(A_mesh-g* R_mesh+c0); dR=makelong*dR
    ax1.streamplot(A_mesh,R_mesh,dA, dR,color='forestgreen')

    camera.snap()    
ani = camera.animate()
   #%%
OutPath = 'C:\\Users\\ellin\\Dropbox\\Prospectus\\Presentation\\nullcline_videos\\'
ani_name = 'Sgro_rap_w_streams'
ani.save(OutPath + ani_name + '.mp4') 
#%% Kamino 2017, vector field with ramp cAMP input
from Kamino2017_agent_and_pop_FUN import Kamino2017_agent

tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
AgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}

dt=0.001;
# initializations
x0=0.01; y0=0.05; z0=0.005
x_trace=[x0]; y_trace=[y0]
Kamino_agent=Kamino2017_agent([x0,y0,z0],AgentParam)

# define extracellular stim trace
signal = 1
Exp_time = Sgro2015Figure3excel["Ramp Input (min Time)"]# np.linspace(0,20,num=1000) # t step 0.02
Exp_time = Exp_time[~np.isnan(Exp_time)]/4*Nt_Kamino
RampInput_Exp = Sgro2015Figure3excel["Ramp Input (nM cAMP)"]
RampInput_Exp = RampInput_Exp[~np.isnan(RampInput_Exp)]
RampInput_Exp = RampInput_Exp[math.floor(len(RampInput_Exp)*4/8): math.floor(len(RampInput_Exp)*7/8)]
Exp_time = Exp_time[0:len(RampInput_Exp)]

Kamino_time = np.arange(0,7.5*Nt_Kamino,dt)
RampInput_Kamino= np.interp(Kamino_time, Exp_time,RampInput_Exp) * signal

for i in range(len(Kamino_time)-1):
    signal_now = RampInput_Kamino[i]
    x_now=x_trace[i]
    y_now=y_trace[i]
    x_next,y_next,z_next= Kamino_agent.update(dt, signal_now)
    x_trace.append(x_next)
    y_trace.append(y_next)

x_trace = np.array(x_trace) 
y_trace = np.array(y_trace)
y_trace_hnorm = (y_trace-Nh_Kamino_offset)/Nh_Kamino
    
t_plot_Kamino = np.array(Kamino_time)/(Nt_Kamino)

# nullclines
x_null_short = np.linspace(-0.5,2,num=100)
dxdt_null_no_stim = 0 + delta
dxdt_null = signal + delta
dydt_null_no_stim_short = (0+delta)**n/((0+delta)**n+(np.power(K*x_null_short,n)))
dydt_null = (signal+delta)**n/((signal+delta)**n+(np.power(K*x_null_short,n)))

#%%
np.savez('figS3_nullcline_ramp_200422.npz',
         t_plot_Goldbeter = t_plot_Goldbeter, 
         b_trace_hnorm = b_trace_hnorm,
         dpdt_null = dpdt_null, dbdt_null = dbdt_null,
         dbdt_null_no_stim = dbdt_null_no_stim,
         p_null = p_null, dpdt_null_no_stim = dpdt_null_no_stim,
         p_trace = p_trace, b_trace = b_trace,
         p_mesh= p_mesh, b_mesh = b_mesh, dp = dp, db = db,
         
         t_plot_Sgro = t_plot_Sgro, A_trace_plot = A_trace_plot,
         A_null = A_null, dRdt_null = dRdt_null,
         dAdt_null = dAdt_null, dAdt_null_no_stim =dAdt_null_no_stim,
         A_trace_orig = A_trace_orig, R_trace_orig =  R_trace_orig,
         A_mesh = A_mesh,R_mesh = R_mesh, dA=dA, dR = dR,
         
         t_plot_Kamino = t_plot_Kamino, y_trace_hnorm = y_trace_hnorm,
         dxdt_null_no_stim = dxdt_null_no_stim, 
         dxdt_null = dxdt_null, dydt_null = dydt_null,
         x_null_short = x_null_short, dydt_null_no_stim_short = dydt_null_no_stim_short,
         x_trace= x_trace, y_trace = y_trace,
         x_mesh = x_mesh, y_mesh = y_mesh, dx = dx, dy=dy)
#%%
#plot out outputs
trace_width = 3
fig5 = plt.figure(figsize=(8,9))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.8)

ax0 = fig5.add_subplot(grid[0, 0],xticks=[0,2.5,5,7.5])
ax0.plot(t_plot_Kamino, y_trace_hnorm, 'green',linewidth=trace_width)
# ax0.axvline(x=5, ls='--', linewidth=trace_width, color='dimgrey')
ax0.axvspan(2.5, 7.5, alpha=0.2, color='b');
ax0.set_title(r'$cAMP_{e}$ ramp input',fontsize=label_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_xlim([0,7.5]); ax0.set_ylim([-0.05,0.4])
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax1= fig5.add_subplot(grid[1:, 0])
ax1.axvline(x=dxdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax1.axvline(x=dxdt_null, ls='-', linewidth=trace_width, color='dimgrey')
ax1.plot(x_null_short,dydt_null_no_stim_short,'lightblue',linewidth=trace_width)
ax1.plot(x_null_short, dydt_null,'deepskyblue',linewidth=trace_width)

ax1.plot( x_trace, y_trace,'green',linewidth=trace_width)
ax1.plot( x_trace_1[0], y_trace_1[0],'*',markersize = 12,color='springgreen')
ax1.set_xlim([-0.2,1.5]); ax1.set_ylim([-0.1,1.1])
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)

# vector field for tiny cAMPe inputs 
x_arr = np.arange(-1,3,0.1)
y_arr = np.arange(-0.5,1.5,0.1)
x_mesh, y_mesh = np.meshgrid(x_arr , y_arr )
makelong = 1
dx = 1/tau*(signal+delta-x_mesh); dx = makelong*dx
dy = (signal+delta)**n/((signal+delta)**n+(np.power(K*x_mesh,n))) - y_mesh;  dy = makelong*dy

ax1.streamplot(x_mesh,y_mesh,dx, dy,color='forestgreen')

fig5.text(0.04, 0.95,'C',ha='center', va='center', fontsize=abcd_font_size)
fig5.text(0.03, 0.8,  r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)
#%%  Make video
fig5 = plt.figure(figsize=(8,7))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.55)
fig5.text(0.04, 0.95,'C',ha='center', va='center', fontsize=abcd_font_size)
fig5.text(0.03, 0.8,  r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.65,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.5, 0.04,'Inhibitor',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.03, 0.35, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

ax0 = fig5.add_subplot(grid[0, 0],xticks=[0,2.5,5,7.5])
ax1= fig5.add_subplot(grid[1:, 0])
camera = Camera(fig5)
down_sample_frame =  500

for huhu in range(int(len(t_plot_Kamino)/down_sample_frame)): 
    i = huhu*down_sample_frame
    ax0.plot(t_plot_Kamino[0:i], y_trace_hnorm[0:i], 'green',linewidth=trace_width)
    ax0.axvspan(2.5, 7.5, alpha=0.2, color='b');
    ax0.set_title(r'$cAMP_{e}$ ramp input',fontsize=label_font_size)
    ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax0.set_xlim([0,7.5]); ax0.set_ylim([-0.05,0.4])
    ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    # nullclines
    signal = RampInput_Kamino[i]
    x_null_short = np.linspace(-0.5,2,num=100)
    dxdt_null = signal + delta
    dydt_null = (signal+delta)**n/((signal+delta)**n+(np.power(K*x_null_short,n)))

    ax1.axvline(x=dxdt_null, ls='-', linewidth=trace_width, color='dimgrey')
    ax1.plot(x_null_short, dydt_null,'deepskyblue',linewidth=trace_width)
    
    ax1.plot( x_trace[0:i], y_trace[0:i],'green',linewidth=trace_width)
    ax1.plot( x_trace[i], y_trace[i],'o',markersize = 8,color='springgreen')
    ax1.set_xlim([-0.2,1.5]); ax1.set_ylim([-0.1,1.1])
    ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    # # vector field for tiny cAMPe inputs 
    # x_arr = np.arange(-1,3,0.1)
    # y_arr = np.arange(-0.5,1.5,0.1)
    # x_mesh, y_mesh = np.meshgrid(x_arr , y_arr )
    # makelong = 1
    # dx = 1/tau*(signal+delta-x_mesh); dx = makelong*dx
    # dy = (signal+delta)**n/((signal+delta)**n+(np.power(K*x_mesh,n))) - y_mesh;  dy = makelong*dy
    # ax1.streamplot(x_mesh,y_mesh,dx, dy,color='forestgreen')
    
    camera.snap()    
ani = camera.animate()
   #%%
OutPath = 'C:\\Users\\ellin\\Dropbox\\Prospectus\\Presentation\\nullcline_videos\\'
ani_name = 'Kamino_rap_wo_streams'
ani.save(OutPath + ani_name + '.mp4') 