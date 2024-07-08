# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 20:45:13 2020

@author: Chuqiao H

PKAR mutant modeling

"""
import os
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


mycwd = os.getcwd()
os.chdir('C:\\Users\\ellin\\Documents\\GitHub\\labscripts\\Plotting Scripts')
from SgroLabPlotStyles import setPlotFont, cm2Inch, sgroStyleDict, greens, blues, pinks, divcolors, oranges
# os.chdir(mycwd)
setPlotFont('Roboto')


#%% WT outputs

from Sgro2015_agent_and_pop_FUN import Sgro2015_agent

from Params import SgroAgentParam
# Sgro 2015 parameters
SgroAgentParam={'e':0.1,'tauA':0.09,'tauR':0.9,'g':0.5,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}

    # PKAR mutant data
a1 = 0.62; b1 = 1.25
SgroAgentParam['g'] = 1/a1
SgroAgentParam['c0'] = b1/a1

A0=-1.5; R0=-0.5
Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)
dt=0.005 ; t_tot=12*Nt_Sgro; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
constant_signal= 10 #np.logspace(-3, 2, num=5)
# constant_signal_arr = np.array([10000])

stim_time_step=int(round(1/12*t_tot/dt)) # at this time step input is applied
signal_trace=np.zeros(len(t))
signal_trace[stim_time_step:] = constant_signal
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
A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
A_trace_plot=(A_trace_orig - Nh_Sgro_offset)/Nh_Sgro;
    
label_fs = 25
trace_width = 6.0 
# colors = plt.cm.summer(np.linspace(0,1,len(constant_signal_arr)+1))
t_plot_Sgro = np.array(t)/(Nt_Sgro)
fig3 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
ax1= fig3.add_subplot(grid[0, 0])

ax1.plot(t_plot_Sgro,A_trace_plot,lw=2.5, color = 'green', label = 'Input='+str(signal))
ax1.axvline(1,color='grey',linewidth='2.5')
ax1.set_ylabel( '\cAMP_{i} trace, A.U.',fontsize=label_fs)
ax1.set_xlabel('Time, A.U.',fontsize=label_fs)
ax1.set_title('Sgro 2015')
leg = ax1.legend();

A_PKAR_10 = A_trace_plot
#%%
np.savez('PKAR_WT_OUT_200809.npz', 
A_WT_1= A_WT_1, A_WT_10k= A_WT_10k,
A_PKAR_10= A_PKAR_10, A_PKAR_10k= A_PKAR_10k)

#%% Load NPZ file
npzname = "PKAR_WT_OUT_200809.npz"
npzfile = np.load( npzname,allow_pickle=True)
A_WT_1 = npzfile['A_WT_1'] ; A_WT_10k = npzfile['A_WT_10k'] 
A_PKAR_1 = npzfile['A_PKAR_1'] ; A_PKAR_10k = npzfile['A_PKAR_10k'] 

#%% Experimental data
# WT
my_dir = '//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/exp_data/'
Sgro2015Figure1excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure1')
#%% Plot results
abcd_font_size = 26
label_font_size=18
title_font_size = 22
sublabel_font_size = 20
trace_width=3
tick_font_size=18

fig5 = plt.figure(figsize=(8, 3.5))
grid = plt.GridSpec(2, 2, wspace=0.2, hspace=0.7)

ax0= fig5.add_subplot(grid[0, 0])
ax0.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 1 FRET Trace (1nM)"],
                              linewidth=trace_width,color='k')
ax0.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 2 FRET Trace (1nM)"],
                               linewidth=trace_width,color='dimgrey')
ax0.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 3 FRET Trace (1nM)"],
                               linewidth=trace_width,color='darkgrey')
#ax1.set_ylabel(r'FRET Signal, A.U.',fontsize=label_font_size)
ax0.set_xlabel('Time (min)',fontsize=label_font_size)
ax0.axvline(x=5, ls='--', linewidth=trace_width, color='green') #dashed line at 5 (cAMP onset)
ax0.set_ylim([-0.1,0.7]); ax0.set_xlim([0, 30])
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_title('1 nM cAMP',color = 'k', fontsize = label_font_size)
ax0.text(-0.26, 0.5, 'Experiment,\nWT FRET ', rotation=90, ha='center',va='center',
         transform = ax0.transAxes, color = 'k', fontsize= tick_font_size)


ax1= fig5.add_subplot(grid[0, 1])
ax1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 1 FRET Trace (10uM)"],
                              linewidth=trace_width,color='k')
ax1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 2 FRET Trace (10uM)"],
                               linewidth=trace_width,color='dimgrey')
ax1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 3 FRET Trace (10uM)"],
                               linewidth=trace_width,color='darkgrey')
ax1.set_xlabel('Time (min)',fontsize=label_font_size)
ax1.axvline(x=5, ls='--', linewidth=trace_width, color='green') #dashed line at 5 (cAMP onset)
ax1.set_ylim([-0.1,0.7]); ax1.set_xlim([0, 30])
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('10 uM cAMP',color = 'k', fontsize = label_font_size)

ax21= fig5.add_subplot(grid[1, 0])
ax21.plot(t_plot_Sgro, A_WT_1, lw = trace_width,color=blues[3])
ax21.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax21.axvline(x=1, ls='--', linewidth=trace_width, color='green') #dashed line at 5 (cAMP onset)
ax21.set_ylim([-0.25,1.1]); ax21.set_xlim([0, 12])
ax21.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax21.text(-0.26, 0.5, 'Model,\nWT ', rotation=90, ha='center',va='center',
         transform = ax21.transAxes, color = 'k', fontsize= tick_font_size)
# ax21.text(-0.28, 1.25, 'B', ha='center',va='center',
#          transform = ax21.transAxes, fontsize=abcd_font_size)

ax22= fig5.add_subplot(grid[1, 1])
ax22.plot(t_plot_Sgro, A_WT_10k, lw=trace_width,color=blues[3])
ax22.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax22.axvline(x=1, ls='--', linewidth=trace_width, color='green') #dashed line at 5 (cAMP onset)
ax22.set_ylim([-0.25,1.1]); ax22.set_xlim([0, 12])
ax22.tick_params(axis='both', which='major', labelsize=tick_font_size)

#%% PKAR mutant exp
abcd_font_size = 26
label_font_size=18
title_font_size = 22
sublabel_font_size = 20
trace_width=3
tick_font_size=18
fig5 = plt.figure(figsize=(8, 3.5))
grid = plt.GridSpec(2, 2, wspace=0.2, hspace=0.7)

t_plot = np.arange(0,160)* 1/3
ax4= fig5.add_subplot(grid[0, 0])
ax4.plot(t_plot, PKAR_exp_10nM[0,:],linewidth=trace_width,color='k')
ax4.plot(t_plot, PKAR_exp_10nM[1,:],linewidth=trace_width,color='dimgrey')
ax4.plot(t_plot, PKAR_exp_10nM[2,:],linewidth=trace_width,color='darkgrey')
ax4.set_xlabel('Time (min)',fontsize=label_font_size)
ax4.axvline(x=5, ls='--', linewidth=trace_width, color='green') #dashed line at 5 (cAMP onset)
ax4.set_ylim([-0.1,0.7]); ax4.set_xlim([0, 40])
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.text(-0.26, 0.5, 'Experiment,\nPKAR FRET ', rotation=90, ha='center',va='center',
         transform = ax4.transAxes, color = 'k', fontsize= tick_font_size)
ax4.set_title('10 nM cAMP',color = 'k', fontsize = label_font_size)

ax41= fig5.add_subplot(grid[0, 1])
ax41.plot(t_plot, PKAR_exp_10uM[0,:],linewidth=trace_width,color='k')
ax41.plot(t_plot, PKAR_exp_10uM[1,:],linewidth=trace_width,color='dimgrey')
ax41.plot(t_plot, PKAR_exp_10uM[2,:],linewidth=trace_width,color='darkgrey')
ax41.set_xlabel('Time (min)',fontsize=label_font_size)
ax41.axvline(x=5, ls='--', linewidth=trace_width, color='green') #dashed line at 5 (cAMP onset)
ax41.set_ylim([-0.1,0.7]); ax41.set_xlim([0, 40])
ax41.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax41.set_title('10 uM cAMP',color = 'k', fontsize = label_font_size)

ax21= fig5.add_subplot(grid[1, 0])
ax21.plot(t_plot_Sgro, A_PKAR_10, lw = trace_width,color=blues[3])
ax21.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax21.axvline(x=1, ls='--', linewidth=trace_width, color='green') #dashed line at 5 (cAMP onset)
ax21.set_ylim([-0.25,1.1]); ax21.set_xlim([0, 12])
ax21.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax21.text(-0.26, 0.5, 'Model,\nPKAR ', rotation=90, ha='center',va='center',
         transform = ax21.transAxes, color = 'k', fontsize= tick_font_size)
# ax21.text(-0.28, 1.25, 'B', ha='center',va='center',
#          transform = ax21.transAxes, fontsize=abcd_font_size)

ax22= fig5.add_subplot(grid[1, 1])
ax22.plot(t_plot_Sgro, A_PKAR_10k, lw=trace_width,color=blues[3])
ax22.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax22.axvline(x=1, ls='--', linewidth=trace_width, color='green') #dashed line at 5 (cAMP onset)
ax22.set_ylim([-0.25,1.1]); ax22.set_xlim([0, 12])
ax22.tick_params(axis='both', which='major', labelsize=tick_font_size)

#%%
import scipy
from scipy.io import loadmat
matdir = 'C:\\Users\\ellin\\Dropbox\\pkaR Mutant Data\\'
annots = loadmat(matdir+'2014-03-28-0X_pkaR-E1c_4uLpmin_10nMcAMPstep_4hr43minStarved.mat')
PKAR_exp_10nM = np.zeros((3,160))
for i in range(3):
    bb = [item.flat[i] for item in annots['outputAll']['output']]
    PKAR_exp_10nM[i,:] = bb[0]
annots = loadmat(matdir+'2014-04-11-60X_pkaR-E1c_4uLpmin_10uMcAMPstep_4hr17minStarved.mat')
PKAR_exp_10uM = np.zeros((5,160))
for i in range(5):
    bb = [item.flat[i] for item in annots['outputAll']['output']]
    PKAR_exp_10uM[i,:] = bb[0]
#%%


#%%
f = h5py.File(matdir+'2014-03-28-0X_pkaR-E1c_4uLpmin_10nMcAMPstep_4hr43minStarved.mat','r')
data = f.get('data/variable1')
data = np.array(data) # For converting to a NumPy array
# %%
PKAR_exp_1 =aa[0,0]