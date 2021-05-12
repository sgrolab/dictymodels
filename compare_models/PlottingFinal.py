# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 22:05:58 2020

@author: cqhuyan

Plotting script for "Bridging Scales to Model Emergent Collective Oscillations in Social Amoeba"

"""
# import packages

import os
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import pandas as pd
import pandas as pd

# Normalization parameters
# from NormParam import *
from NB_pop_functions import plot_POP_oscillation
#Normalization parameters
from Params import NormParams
for key,val in NormParams.items():
        exec(key + '=val')
        
from NB_pop_functions import * 
from iocustom import import_npz

# set up new default font
import matplotlib
font = {'family' : 'Arial'}
matplotlib.rc('font', **font)

def cm2inch(value):
    return value/2.54

#%% imaging colors
colors = np.array([[240,145,147],
                  [240,206,109],
                  [130,191,91],
                  [150,110,184],
                  [125,208,227],
                  [59,34,255], # color for experiments
                  [35,145,40]])# color for models
goldcolor = (colors[0,:]+(255-colors[0,:])*0)/255 # tint factor, larger means more "whiteness" to the original color
maedacolor = (colors[1,:]+(255-colors[1,:])*0)/255
gregorcolor = (colors[2,:]+(255-colors[2,:])*0)/255
sgrocolor = (colors[3,:]+(255-colors[3,:])*0)/255
kaminocolor = (colors[4,:]+(255-colors[4,:])*0)/255

# color combo 2
# colors = np.array([[168,69,154],
#                   [129,129,129],
#                   [87,96,171],
#                   [90,173,90],
#                   [105,201,205],
#                   [59,34,255], # color for experiments
#                   [35,145,40]])# color for models
# goldcolor = (colors[0,:]+(255-colors[0,:])*0)/255 # tint factor, larger means more "whiteness" to the original color
# maedacolor = (colors[1,:]+(255-colors[1,:])*0)/255
# gregorcolor = (colors[2,:]+(255-colors[2,:])*0)/255
# sgrocolor = (colors[3,:]+(255-colors[3,:])*0)/255
# kaminocolor = (colors[4,:]+(255-colors[4,:])*0)/255

expcolor = (colors[5,:]+(255-colors[5,:])*0)/255
simcolor = (colors[6,:]+(255-colors[6,:])*0)/255

#%% Figure 2- single cell adaptive spiking and oscillations
# experimental data
my_dir = r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/'
Sgro2015Figure1excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure1')

# load saved npz output file- adaptive spiking
import_npz('../model_outputs/single_cell_adaptive_200319.npz',globals())

# load saved npz output file- oscillations
npzfile = np.load('../model_outputs/single_cell_oscillation_200319.npz')
# import_npz('single_cell_oscillation_200319.npz',globals())
t_plot_Goldbeter_osc =  npzfile['t_plot_Goldbeter'] ; b_trace_osc = npzfile['b_trace']
t_plot_Maeda_osc=npzfile['t_plot_Maeda'] ; cAMPi_trace_osc=npzfile['cAMPi_trace']
t_plot_Gregor_osc=npzfile['t_plot_Gregor']; gregor_campCyto_trace_osc=npzfile['gregor_campCyto_trace']
t_plot_Sgro_osc = npzfile['t_plot_Sgro']; A_trace_plot_osc = npzfile['A_trace_plot']
t_plot_Kamino_osc = npzfile['t_plot_Kamino']; y_trace_osc = npzfile['y_trace']

# plotting font sizes
abcd_font_size = 28
label_font_size=22
title_font_size = 24
sublabel_font_size = 20
trace_width=3
tick_font_size=22

fig5 = plt.figure(figsize=(cm2inch(45),cm2inch(30)))
grid = plt.GridSpec(9, 2)
# grid = plt.GridSpec(9, 2, wspace=0.3, hspace=1.8)

ax01= fig5.add_subplot(grid[0:2, 0])
ax01.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 1 FRET Trace (1nM)"],
                              linewidth=trace_width,color='k')
ax01.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 2 FRET Trace (1nM)"],
                               linewidth=trace_width,color='dimgrey')
ax01.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 3 FRET Trace (1nM)"],
                               linewidth=trace_width,color='darkgrey')
# ax01.set_ylabel(r'FRET Signal, A.U.',fontsize=label_font_size)
ax01.set_xlabel('Time (min)',fontsize=label_font_size)
ax01.axvline(x=5, ls='--', linewidth=trace_width, color = expcolor) #dashed line at 5 (cAMP onset)
ax01.set_ylim([-0.1,0.7]); ax01.set_xlim([0, 30])
ax01.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax01.set_title('Experiment',color = 'k', fontsize = title_font_size)
ax01.text(0.7,0.8,' 1nM cAMP', horizontalalignment='center',verticalalignment='center',
         transform = ax01.transAxes, color = 'k', fontsize=label_font_size)
ax01.text(-0.19, 1.25, 'A',
         horizontalalignment='center',verticalalignment='center',
         transform = ax01.transAxes, color = expcolor, fontsize=abcd_font_size)

ax02= fig5.add_subplot(grid[3:7, 0])
#ax2.axvline(x = 5 , ymin=-0.2, ymax = 1.2, ls = '--', 
#            linewidth=trace_width, color=mycolors[6])
ax02.plot(t_plot_Goldbeter, b_trace,linewidth=trace_width, color = goldcolor,
         label='Receptor desensitization')
ax02.plot(t_plot_Maeda, cAMPi_trace,linewidth=trace_width,color = maedacolor,
         label='CDINFB')
ax02.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax02.plot(t_plot_Gregor,gregor_campCyto_trace,linewidth=trace_width, 
         color= gregorcolor,label='Phase oscillator')
ax02.plot(t_plot_Sgro, A_trace_plot,linewidth=trace_width,
         color=sgrocolor, label='IPNFB')
ax02.plot(t_plot_Kamino, y_trace,linewidth=trace_width, 
         ls = '--',color=kaminocolor, label='IFFL')
ax02.set_title('Simulation',color = 'k', fontsize = title_font_size)
ax02.text(-0.19, 1.05, 'C',
         horizontalalignment='center',verticalalignment='center',
         transform = ax02.transAxes, color = simcolor, fontsize=abcd_font_size)
ax02.set_ylim([-0.4,1.2])#([-0.2,1.2]) # adaptive spikes
ax02.set_xlim([0,6])
# ax02.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax02.tick_params(axis='both', which='major', labelsize=tick_font_size)

# leg = ax2.legend()
ax02.legend(bbox_to_anchor=(2.1,-0.28),ncol=5,prop={'size': 18})


ax1= fig5.add_subplot(grid[0:2, 1])
ax1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 1 FRET Trace (10uM)"],
                              linewidth=trace_width,color='k')
ax1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 2 FRET Trace (10uM)"],
                               linewidth=trace_width,color='dimgrey')
ax1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 3 FRET Trace (10uM)"],
                               linewidth=trace_width,color='darkgrey')

ax1.set_xlabel('Time (min)',fontsize=label_font_size)
ax1.axvline(x=5, ls='--', linewidth=trace_width, color = expcolor) #dashed line at 5 (cAMP onset)
ax1.set_ylim([-0.1,0.7]); ax1.set_xlim([0, 30])
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Experiment',color = 'k', fontsize = title_font_size)
ax1.text(0.7,0.8,' 10uM cAMP', horizontalalignment='center',verticalalignment='center',
         transform = ax1.transAxes, color = 'k', fontsize=label_font_size)
ax1.text(-0.13, 1.25, 'B',
         horizontalalignment='center',verticalalignment='center',
         transform = ax1.transAxes, color = expcolor, fontsize=abcd_font_size)

ax2= fig5.add_subplot(grid[3:7, 1])
#ax2.axvline(x = 5 , ymin=-0.2, ymax = 1.2, ls = '--', 
#            linewidth=trace_width, color=mycolors[6])
ax2.plot(t_plot_Goldbeter_osc, b_trace_osc,linewidth=trace_width, color=goldcolor,
         label='Receptor desensitization')
ax3 = ax2.twinx()
ax3.plot(t_plot_Maeda_osc, cAMPi_trace_osc,linewidth=trace_width,color=maedacolor,
         label='CDINFB')
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_ylabel(r'CDINFB $cAMP_{i}$, A.U.', color=maedacolor,
                fontsize=label_font_size-3)
ax3.set_ylim([-55,260])

ax2.plot(t_plot_Gregor_osc,gregor_campCyto_trace_osc,linewidth=trace_width, 
         color = gregorcolor,label='Phase oscillator')
ax2.plot(t_plot_Sgro_osc, A_trace_plot_osc,linewidth=trace_width,
         color= sgrocolor, label='IPNFB')
ax2.plot(t_plot_Kamino_osc, y_trace_osc, linewidth=trace_width, 
         ls = '--',color=kaminocolor, label='IFFL')
ax2.set_title('Simulation',color = 'k', fontsize = title_font_size)
ax2.text(-0.13, 1.05, 'D',
         ha='center',va='center',transform = ax2.transAxes, 
         color = simcolor, fontsize=abcd_font_size)

ax2.set_ylim([-0.4,1.85])# oscillations
ax2.set_xlim([0,6])
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)

fig5.text(0.08, 0.82, 'FRET, A.U.', ha='center', va='center',rotation='vertical',
          fontsize=label_font_size)
fig5.text(0.08, 0.45,r'$cAMP_{i}$, A.U.', ha='center', va='center',rotation='vertical',
          fontsize=label_font_size)
fig5.text(0.5, 0.2, 'Time, A.U.', ha='center', va='center',fontsize=label_font_size)

grid.tight_layout(fig5,rect=[0, 0, 1, 1],pad = 10)
# fig5.tight_layout()
plt.show()

fig5.savefig('Fig2_20200530_orig.png', bbox_inches='tight')#
#%% Figure S1- nullclines for bifurcation analysis
# import simulation outputs
import_npz('../model_outputs/figS1_nullcline_200422.npz',globals())

abcd_font_size = 28
label_font_size=20
title_font_size = 22
sublabel_font_size = 18
trace_width=3
tick_font_size=20
m_size = 18 # marker size
m_color = np.array([0,255,127])/255

fig5 = plt.figure(figsize=(cm2inch(40),cm2inch(32)))
grid = plt.GridSpec(7, 4, wspace=0.4, hspace=1)
# Goldbeter 1987
ax0 = fig5.add_subplot(grid[0, 0])
ax0.plot(t_plot_Goldbeter, b_trace_1_hnorm, color = 'darkgreen',linewidth=trace_width)
ax0.axvline(x=1, ls='--', linewidth=trace_width, color='darkgrey')
ax0.set_title(r'$cAMP_{e}$ input: '+str(signal_1),fontsize=sublabel_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_xlim([0,6]); ax0.set_ylim([-0.2,2])
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_ylabel(r'$cAMP_{i}$', fontsize=sublabel_font_size)


ax1= fig5.add_subplot(grid[1:3, 0])
ax1.axvline(x=dpdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax1.axvline(x=dpdt_null_1, ls='-', linewidth=trace_width, color='dimgrey')
ax1.plot(p_null,dbdt_null_no_stim, 'lightblue',linewidth=trace_width)
ax1.plot(p_null, dbdt_null_1,'deepskyblue',linewidth=trace_width)
ax1.plot( p_trace_1, b_trace_1,color = 'darkgreen',linewidth= trace_width)
ax1.plot(p_trace_1[0],b_trace_1[0],'*',markersize = m_size,color= m_color)
ax1.set_xlim([-0.5,1.5]); ax1.set_ylim([-20,700])
# ax1.set_ylim([0,0.5]); 
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_ylabel(r'$cAMP_{i}$', fontsize=sublabel_font_size)

ax00 = fig5.add_subplot(grid[0, 1])
ax00.plot(t_plot_Goldbeter, b_trace_10k_hnorm, color = 'darkgreen',linewidth=trace_width)
ax00.axvline(x=1, ls='--', linewidth=trace_width, color = 'darkgrey')
ax00.set_title(r'$cAMP_{e}$ input: '+str(signal_10k),fontsize= sublabel_font_size)
ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax00.set_xlim([0,6]); ax00.set_ylim([-0.2,2])
ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax2= fig5.add_subplot(grid[1:3, 1])
ax2.axvline(x=dpdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax2.axvline(x=dpdt_null_10k, ls='-', linewidth=trace_width, color='dimgrey')
ax2.plot(p_null,dbdt_null_no_stim, 'lightblue',linewidth=trace_width)
ax2.plot(p_null, dbdt_null_10k,'deepskyblue',linewidth=trace_width)
ax2.plot( p_trace_10k, b_trace_10k,color = 'darkgreen',linewidth= 4)
ax2.plot(p_trace_10k[0],b_trace_10k[0],'*',markersize = m_size,color=m_color)
ax2.set_xlim([-0.5,1.5]); ax2.set_ylim([-20,700])
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax0.text(-0.4, 2.1, 'A',
         ha='center',va='center',transform = ax0.transAxes, fontsize=abcd_font_size)
ax0.text(1.16, -0.8,'Time, A.U.',ha='center', va='center',
         transform = ax0.transAxes, fontsize=sublabel_font_size)
ax1.text(1.16, -0.3,r'Inhibitor',ha='center', va='center',
         transform = ax1.transAxes, fontsize=sublabel_font_size)
ax0.text(1.16, 2,'Receptor desensitization',ha='center', va='center',
         transform = ax0.transAxes,color= goldcolor,fontsize=title_font_size)

# IPNFB
ax3 = fig5.add_subplot(grid[0, 2])
ax3.plot(t_plot_Sgro, A_trace_plot_1, color = sgrocolor ,linewidth=trace_width)
ax3.axvline(x=1, ls='--', linewidth=trace_width, color='darkgrey')
ax3.set_title(r'$cAMP_{e}$ input: '+str(signal_1),fontsize=sublabel_font_size)
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_xlim([0,6]); ax3.set_ylim([-0.4,1.2])
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_ylabel(r'$cAMP_{i}$'+'\n(Activator)', fontsize=sublabel_font_size)


ax4= fig5.add_subplot(grid[1:3, 2])
ax4.plot(A_null,dRdt_null, 'darkgrey',linewidth=trace_width)
ax4.plot(A_null,dAdt_null_no_stim,'lightblue',linewidth=trace_width)
ax4.plot(A_null, dAdt_null_1,'deepskyblue',linewidth=trace_width)
ax4.plot( A_trace_orig_1, R_trace_orig_1, color = sgrocolor ,linewidth=trace_width)
ax4.plot(A_trace_orig_1[0],R_trace_orig_1[0],'*',markersize = m_size,color=m_color)
ax4.set_ylim([-1,2.5])
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_ylabel('Inhibitor', fontsize=sublabel_font_size)

ax5 = fig5.add_subplot(grid[0, 3])
ax5.plot(t_plot_Sgro, A_trace_plot_10k, color = sgrocolor ,linewidth=trace_width)
ax5.axvline(x=1, ls='--', linewidth=trace_width, color= 'darkgrey')
ax5.set_title(r'$cAMP_{e}$ input: '+str(signal_10k),fontsize = sublabel_font_size)
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5.set_xlim([0,6]); ax5.set_ylim([-0.4,1.2])
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax6= fig5.add_subplot(grid[1:3, 3])
ax6.plot(A_null,dRdt_null, 'darkgrey',linewidth=trace_width)
ax6.plot(A_null,dAdt_null_no_stim,'lightblue',linewidth=trace_width)
ax6.plot(A_null, dAdt_null_10k,'deepskyblue',linewidth=trace_width)
ax6.plot( A_trace_orig_10k, R_trace_orig_10k, color = sgrocolor ,linewidth=trace_width)
ax6.plot(A_trace_orig_10k[0],R_trace_orig_10k[0],'*',markersize = m_size,color=m_color)
ax6.set_ylim([-1,2.5])
ax6.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax3.text(-0.4, 2.1, 'B',
         ha='center',va='center',transform = ax3.transAxes, fontsize=abcd_font_size)
ax3.text(1.16, -0.8,'Time, A.U.',ha='center', va='center',
         transform = ax3.transAxes, fontsize=sublabel_font_size)
ax4.text(1.16, -0.35,r'$cAMP_{i}$'+'(Activator)',ha='center', va='center',
         transform = ax4.transAxes, fontsize=sublabel_font_size)
ax3.text(1.16, 2,'IPNFB',ha='center', va='center',
         transform = ax3.transAxes,color= sgrocolor,fontsize=title_font_size)

# IFFL
ax7 = fig5.add_subplot(grid[4, 0])
ax7.plot(t_plot_Kamino, y_trace_1_hnorm, color = kaminocolor,linewidth=trace_width)
ax7.axvline(x=1, ls='--', linewidth=trace_width, color='darkgrey')
ax7.set_title(r'$cAMP_{e}$ input: '+str(signal_1),fontsize=label_font_size)
ax7.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax7.set_xlim([0,6]); ax7.set_ylim([-0.2,1.2])
ax7.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax7.set_ylabel(r'$cAMP_{i}$', fontsize=sublabel_font_size)


ax8= fig5.add_subplot(grid[5:, 0])
ax8.axvline(x=dxdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax8.axvline(x=dxdt_null_1, ls='-', linewidth=trace_width, color='dimgrey')
ax8.plot(x_null_short,dydt_null_no_stim_short,'lightblue',linewidth=trace_width)
ax8.plot(x_null_short, dydt_null_1,'deepskyblue',linewidth=trace_width)

ax8.plot( x_trace_1, y_trace_1,color = kaminocolor, linewidth=trace_width)
ax8.plot( x_trace_1[0], y_trace_1[0],'*',markersize = m_size,color=m_color)
ax8.set_xlim([-0.5,2])
ax8.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax8.set_ylabel(r'$cAMP_{i}$', fontsize=sublabel_font_size)

ax9 = fig5.add_subplot(grid[4, 1])
ax9.plot(t_plot_Kamino, y_trace_10k_hnorm, color = kaminocolor, linewidth=trace_width)
ax9.axvline(x=1, ls='--', linewidth=trace_width, color='darkgrey')
ax9.set_title(r'$cAMP_{e}$ input: '+str(signal_10k),fontsize=label_font_size)
ax9.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax9.set_xlim([0,6]); ax9.set_ylim([-0.2,1.2])
ax9.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax10= fig5.add_subplot(grid[5:, 1])
ax10.axvline(x=dxdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax10.axvline(x=dxdt_null_10k, ls='-', linewidth=trace_width, color='dimgrey')
ax10.plot(x_null_long,dydt_null_no_stim_long,'lightblue',linewidth=trace_width)
ax10.plot(x_null_long, dydt_null_10k,'deepskyblue',linewidth=trace_width)
ax10.plot( x_trace_10k, y_trace_10k,color = kaminocolor,linewidth=trace_width)
ax10.plot( x_trace_10k[0], y_trace_10k[0],'*',markersize = m_size,color=m_color)
ax10.set_xlim([-1000,11000])
ax10.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax7.text(-0.4, 2.1, 'C',
         ha='center',va='center',transform = ax7.transAxes, fontsize=abcd_font_size)
ax7.text(1.16, -0.8,'Time, A.U.',ha='center', va='center',
         transform = ax7.transAxes, fontsize=sublabel_font_size)
ax8.text(1.16, -0.3,r'Inhibitor',ha='center', va='center',
         transform = ax8.transAxes, fontsize=sublabel_font_size)
ax7.text(1.16, 2,'IFFL',ha='center', va='center',
         transform = ax7.transAxes,color= kaminocolor,fontsize=title_font_size)

grid.tight_layout(fig5,rect=[0, 0, 1, 1],pad = 10)

fig5.savefig('figS1_nullclines_tight_200530.png', bbox_inches='tight')
#%% Figure S2- excitability
import_npz('../model_outputs/figS2_excitability_200422.npz',globals())

m_size = 20 # marker size
m_color = np.array([0,255,127])/255

fig5 = plt.figure(figsize=(cm2inch(24),cm2inch(18)))
grid = plt.GridSpec(3, 2, wspace=0.3, hspace=0.8)

ax0 = fig5.add_subplot(grid[0, 0])
ax0.plot(t_plot_Sgro, A_trace_plot_smalle, color = sgrocolor, linewidth=trace_width)
ax0.axvline(x=1, ls='--', linewidth=trace_width, color='k')
ax0.set_title('Excitability: '+str(e_small),fontsize=label_font_size)
ax0.set_xlim([0,6]); ax0.set_ylim([-0.4,1.2])
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax1= fig5.add_subplot(grid[1:, 0])
ax1.plot( A_null,dRdt_null,'darkgrey',linewidth=trace_width)
ax1.plot( A_null,dAdt_null,'deepskyblue',linewidth=trace_width)
ax1.plot( A_trace_orig_smalle, R_trace_orig_smalle, color = sgrocolor,linewidth=trace_width)
ax1.set_ylim([-0.7,3])
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax1.streamplot(A_mesh,R_mesh,dA, dR,density=1,color='forestgreen')
ax1.plot(A_trace_orig_smalle[0],R_trace_orig_smalle[0],'*',markersize =  m_size,color = m_color)

ax00 = fig5.add_subplot(grid[0, 1])
ax00.plot(t_plot_Sgro, A_trace_plot_bige,  color = sgrocolor, linewidth=trace_width)
ax00.axvline(x=1, ls='--', linewidth=trace_width, color= 'k')
ax00.set_title('Excitability: '+str(e_big),fontsize=label_font_size)
ax00.set_xlim([0,6]); ax00.set_ylim([-0.4,1.2])
ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax2= fig5.add_subplot(grid[1:, 1])
ax2.plot(A_null,dRdt_null, 'darkgrey',linewidth=trace_width)
ax2.plot(A_null, dAdt_null,'deepskyblue',linewidth=trace_width)
ax2.plot( A_trace_orig_bige, R_trace_orig_bige, color = sgrocolor, linewidth=trace_width)
ax2.set_ylim([-0.7,3])
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax2.streamplot(A_mesh, R_mesh,dA,dR, density=1,color='forestgreen')
ax2.plot(A_trace_orig_bige[0],R_trace_orig_bige[0],'*',markersize = m_size,color= m_color)

fig5.text(0.04, 0.95,'A',ha='center', va='center', fontsize=abcd_font_size)
fig5.text(0.04, 0.8, r'$cAMP_{i}$'+'\n(Activator)', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
fig5.text(0.5, 0.62,'Time, A.U.',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.5, 0.02,r'$cAMP_{i}$'+'(Activator)',ha='center', va='center',fontsize=label_font_size)
fig5.text(0.04, 0.35, 'Inhibitor', ha='center', va='center', rotation='vertical',fontsize= label_font_size)

fig5.text(0.5, 0.95,'IPNFB',ha='center', va='center',color=sgrocolor,fontsize=title_font_size)

grid.tight_layout(fig5,rect=[0, 0, 1, 1],pad = 10) # can not make all axes height small enough to accomodate all axes decorations 
plt.show()

fig5.savefig('figS2_Sgro_excitability_tight_200530.png', bbox_inches='tight')

#%% Fig 3 single cell step ramp
# experimental data
my_dir = r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/'
Sgro2015Figure3excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure3')
# load simulation data
import_npz('../model_outputs/single_cell_StepRamp_200320.npz',globals())

abcd_font_size = 28
label_font_size=20
title_font_size = 22
sublabel_font_size = 18
trace_width=3
tick_font_size=20

fig3 = plt.figure(figsize=(cm2inch(35),cm2inch(18)))
grid = plt.GridSpec(4, 3, wspace=0.31, hspace=0.9)

# ax0 = fig3.add_subplot(grid[0, 0])
ax0 = fig3.add_axes([0.085, 0.72, 0.245, 0.1])
ax0.plot(Sgro2015Figure3excel["Ramp Input (min Time)"],Sgro2015Figure3excel["Ramp Input (nM cAMP)"],
                              color='k', linewidth = trace_width)
ax0.set_title('Experiment',color='k', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax0.set_ylabel(r'$cAMP_{e}$'+'\n(nM)',fontsize=sublabel_font_size)
#ax0.text(-16, 0.5, r'$cAMP_{e}$(nM)', ha='center',va='center',rotation ='vertical',
#         color = 'k', fontsize=label_font_size)
# ax0.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_xlim([0,80])
ax0.set_ylim([-0.1,1.1])
ax0.xaxis.set_major_formatter(plt.NullFormatter()) # hide x axis
ax0.text(-0.33, 1.6, 'A', ha='center',va='center',
         transform = ax0.transAxes, color = 'k', fontsize=abcd_font_size)

ax1= fig3.add_subplot(grid[1:3, 0])
ax1.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 1 FRET Trace"],
                               color='k', linewidth = trace_width)
ax1.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 2 FRET Trace"],
                               color='grey', linewidth = trace_width)
ax1.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 3 FRET Trace"],
                               color='lightgrey', linewidth = trace_width)
ax1.axvspan(10, 30, alpha=0.2, color=expcolor); ax1.axvspan(50, 70, alpha=0.2, color='b')
ax1.set_ylabel('FRET,\n A.U.',fontsize=sublabel_font_size)
ax1.set_xlabel('Time (min)',fontsize=sublabel_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_ylim([-0.1,0.6]); ax1.set_xlim([0,80])

ax2= fig3.add_subplot(grid[0:2, 1])
ax2.plot(t_plot_Goldbeter, b_trace, color=goldcolor,linewidth=trace_width)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
# ax2.set_title('Receptor desensitization',color = goldcolor,  fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax2.set_title('Receptor Desensitization',color = goldcolor,  fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax2.axvspan(2, 6, alpha=0.3, color=simcolor); ax2.axvspan(10, 14, alpha=0.3, color='g')
ax2.set_ylim([-0.2,1.2]); ax2.set_xlim([0,16])
ax2.text(-0.3, 1.14, 'B', ha='center',va='center',
         transform = ax2.transAxes, color = 'k', fontsize=abcd_font_size)

ax3= fig3.add_subplot(grid[0:2, 2])
ax3.plot(t_plot_Maeda, cAMPi_trace, color= maedacolor,linewidth=trace_width)
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
# ax3.set_title('CDINFB', color = maedacolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax3.set_title('CDINFB', color = maedacolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax3.axvspan(2,6, alpha=0.3, color=simcolor); ax3.axvspan(10,14, alpha=0.3, color='g')
ax3.set_ylim([-0.2,1.2]);  ax3.set_xlim([0,16])
# ax3.text(-0.12, 1.06, 'C', ha='center',va='center',
#          transform = ax3.transAxes, color = simcolor, fontsize=abcd_font_size)

ax4= fig3.add_subplot(grid[2:, 1])
ax4.plot(t_plot_Sgro, A_trace_plot, color=sgrocolor,linewidth=trace_width)
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('IPNFB', color = sgrocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax4.axvspan(2,6, alpha=0.3, color=simcolor); ax4.axvspan(10,14, alpha=0.3, color='g')
ax4.set_ylim([-0.4,1.2]);  ax4.set_xlim([0,16])
# ax4.text(-0.12, 1.06, 'D', ha='center',va='center',
#          transform = ax4.transAxes, color = simcolor, fontsize=abcd_font_size)

ax5= fig3.add_subplot(grid[2:, 2])
ax5.plot(t_plot_Kamino, y_trace, color=kaminocolor,linewidth=trace_width)
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5.set_title('IFFL', color = kaminocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax5.axvspan(2,6, alpha=0.3, color= simcolor); ax5.axvspan(10,14, alpha=0.3, color='g')
ax5.set_ylim([-0.4,1.2]);  ax5.set_xlim([0,16])
# ax5.text(-0.12, 1.06, 'E', ha='center',va='center',
#          transform = ax5.transAxes, color = simcolor, fontsize=abcd_font_size)

fig3.text(0.36, 0.5, r'$cAMP_{i}$', fontsize=label_font_size,va='center', rotation='vertical')
fig3.text(0.65, 0.02, 'Time, A.U.', fontsize=label_font_size, ha='center')
grid.tight_layout(fig3,rect=[0, 0, 1, 1]) # can not make all axes height small enough to accomodate all axes decorations 
# plt.show()
# fig3.savefig('fig3_sc_step_ramp_tight_200530.png', bbox_inches='tight')

#%% fig S3 single cell  ramp input nullcline
# load simulation data
import_npz('../model_outputs/figS3_nullcline_ramp_200422.npz',globals())

abcd_font_size = 28
label_font_size=22
title_font_size = 24
sublabel_font_size = 18
trace_width=3
tick_font_size=20
m_size = 18 # marker size
m_color = np.array([0,255,127])/255

fig5 = plt.figure(figsize=(cm2inch(45),cm2inch(20)))
grid = plt.GridSpec(3, 3, wspace=0.3, hspace=0.6)

ax0 = fig5.add_subplot(grid[0, 0],xticks=[0,2.5,5,7.5])
ax0.plot(t_plot_Goldbeter, b_trace_hnorm, color = 'darkgreen',linewidth=trace_width)
ax0.axvspan(2.5, 7.5, alpha=0.2, color= simcolor);
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_xlim([0,7.5]);  ax0.set_ylim([-0.03,0.22])
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax0.set_ylabel(r'$cAMP_{i}$', fontsize= label_font_size)
# ax0.set_xlabel('Time, A.U.',fontsize= label_font_size)
ax0.text(0.68, 0.8, r'$cAMP_{e}$ ramp input', ha='center',va='center',
         transform = ax0.transAxes, fontsize=sublabel_font_size)

ax1= fig5.add_subplot(grid[1:, 0])
ax1.axvline(x=dpdt_null, ls='-', linewidth=trace_width, color='dimgrey')
ax1.plot(p_null, dbdt_null,'deepskyblue',linewidth=trace_width)
ax1.axvline(x=dpdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax1.axvline(x=dpdt_null, ls='-', linewidth=trace_width, color='dimgrey')
ax1.plot(p_null,dbdt_null_no_stim, 'lightblue',linewidth=trace_width)
ax1.plot(p_null, dbdt_null,'deepskyblue',linewidth=trace_width)

ax1.plot( p_trace, b_trace, color = 'darkgreen',linewidth= trace_width)
ax1.plot(p_trace[0],b_trace[0],'*',markersize = m_size,color=m_color)
ax1.set_xlim([-0.5,1.5]); # ax1.set_ylim([-0.5,2])
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax1.set_ylabel(r'$cAMP_{i}$', fontsize= label_font_size)
ax1.set_xlabel('Inhibitor',fontsize= label_font_size)
# vector field for cAMPe = 1
ax1.streamplot(p_mesh,b_mesh,dp, db,color='forestgreen')

ax0.text(-0.24, 1.5, 'A', ha='center',va='center',
         transform = ax0.transAxes, color = 'k', fontsize=abcd_font_size)
ax0.text(0.5, 1.18,'Receptor desensitization',ha='center', va='center',
         transform = ax0.transAxes, color= goldcolor, fontsize=title_font_size)

# IPNFB
ax2 = fig5.add_subplot(grid[0, 1],xticks=[0,2.5,5,7.5])
ax2.plot(t_plot_Sgro, A_trace_plot, color = sgrocolor,linewidth=trace_width)
ax2.axvspan(2.5, 7.5, alpha=0.2, color=simcolor);
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_xlim([0,7.5]); ax2.set_ylim([-0.2,1.2])
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.text(0.68, 0.8, r'$cAMP_{e}$ ramp input', ha='center',va='center',
         transform = ax2.transAxes, fontsize=sublabel_font_size)
ax2.set_xlabel('Time, A.U.',fontsize= label_font_size)

ax3= fig5.add_subplot(grid[1:, 1])
ax3.plot(A_null, dRdt_null,'dimgrey',linewidth=trace_width)
ax3.plot(A_null, dAdt_null,'deepskyblue',linewidth=trace_width)
ax3.plot(A_null,dAdt_null_no_stim,'lightblue',linewidth=trace_width)
ax3.plot( A_trace_orig, R_trace_orig, color = sgrocolor,linewidth= 4)
ax3.plot(A_trace_orig[0],R_trace_orig[0],'*',markersize = m_size,color= m_color)
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_xlim([-2.5,2.5]); ax3.set_ylim([-1,3])
# vector field for tiny cAMPe inputs 
ax3.streamplot(A_mesh,R_mesh,dA, dR,color='forestgreen')
ax3.set_xlabel(r'$cAMP_{i}$ (activator)', fontsize= label_font_size)
ax3.set_ylabel('Inhibitor',fontsize= label_font_size)

ax2.text(-0.24, 1.5, 'B', ha='center',va='center',
         transform = ax2.transAxes, color = 'k', fontsize=abcd_font_size)
ax2.text(0.5, 1.18,'IPNFB',ha='center', va='center',
         transform = ax2.transAxes, color= sgrocolor, fontsize=title_font_size)

# IFFL
ax4= fig5.add_subplot(grid[0, 2],xticks=[0,2.5,5,7.5])
ax4.plot(t_plot_Kamino, y_trace_hnorm, color = kaminocolor,linewidth=trace_width)
ax4.axvspan(2.5, 7.5, alpha=0.2, color=simcolor);
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_xlim([0,7.5]); ax4.set_ylim([-0.05,0.5])
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.text(0.68, 0.8, r'$cAMP_{e}$ ramp input', ha='center',va='center',
         transform = ax4.transAxes, fontsize=sublabel_font_size)

ax5= fig5.add_subplot(grid[1:, 2])
ax5.axvline(x=dxdt_null_no_stim, ls='-', linewidth=trace_width, color='darkgrey')
ax5.axvline(x=dxdt_null, ls='-', linewidth=trace_width, color='dimgrey')
ax5.plot(x_null_short,dydt_null_no_stim_short,'lightblue',linewidth=trace_width)
ax5.plot(x_null_short, dydt_null,'deepskyblue',linewidth=trace_width)

ax5.plot( x_trace, y_trace,color = kaminocolor,linewidth=trace_width)
ax5.plot( x_trace[0], y_trace[0],'*',markersize =m_size,color=m_color)
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5.set_ylabel(r'$cAMP_{i}$', fontsize= label_font_size)
ax5.set_xlabel('Inhibitor',fontsize= label_font_size)
# vector field for tiny cAMPe inputs 
ax5.streamplot(x_mesh,y_mesh,dx, dy,color='forestgreen')
ax5.set_xlim([-0.5,2]); ax5.set_ylim([-0.1,1.1])
ax4.text(-0.24, 1.5, 'C', ha='center',va='center',
         transform = ax4.transAxes, color = 'k', fontsize=abcd_font_size)
ax4.text(0.5, 1.18,'IFFL',ha='center', va='center',
         transform = ax4.transAxes, color= kaminocolor, fontsize=title_font_size)

fig5.savefig('figS3_ramp_nullclines_tight_200530.png', bbox_inches='tight')

#%% Fig 4 single cell entrainment quality
#experimental data
Sgro2015Figure4 = scipy.io.loadmat('C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/Figure4Data.mat')
Sgro2015Figure4excel = pd.read_excel(r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure4')

# load simulation data
import_npz('../model_outputs/single_cell_entrainment_200320.npz',globals())

abcd_font_size = 28
label_font_size=22
title_font_size = 20
sublabel_font_size = 18
trace_width=3
tick_font_size=20

fig3 = plt.figure(figsize=(cm2inch(22),cm2inch(33)))
grid = plt.GridSpec(6, 2, wspace=0.35, hspace=1.5)

ax0u= fig3.add_subplot(grid[0, 0])
ax0u.plot(Sgro2015Figure4excel["Time (min)"]-1,Sgro2015Figure4excel["Cell 1 FRET Trace (1 min pulse)"],
                               'k',linewidth=trace_width)
ax0u.plot(Sgro2015Figure4excel["Time (min)"]-1,Sgro2015Figure4excel["Cell 2 FRET Trace (1 min pulse)"],
                               'dimgrey',linewidth=trace_width)
ax0u.plot(Sgro2015Figure4excel["Time (min)"]-1,Sgro2015Figure4excel["Cell 3 FRET Trace (1 min pulse)"],
                               'darkgrey',linewidth=trace_width)
ax0u.set_xlim([0,23]); ax0u.tick_params(axis='both', which='major', labelsize=tick_font_size)
for i in range(4):
    ax0u.axvspan(4.5+i*6, 4.5+i*6+1, alpha=0.2, color=expcolor)
ax0u.text(-0.35 , 1.85, 'A', ha='center',va='center',
     transform = ax0u.transAxes, color = 'k', fontsize=abcd_font_size)
ax0u.text(-0.3, -0.5, 'FRET, A.U.',ha='center',va='center', rotation='vertical',
     transform = ax0u.transAxes, color = 'k', fontsize=sublabel_font_size)
ax0u.text(0.5, 1.23, r'$cAMP_{e}$ input',ha='center',va='center', 
     transform = ax0u.transAxes, alpha=0.6, color = expcolor, fontsize= sublabel_font_size)

#ax0l= fig3.add_subplot(grid[1, 0])
ax0l = fig3.add_axes([0.13, 0.7, 0.325,0.062])
ax0l.plot(Sgro2015Figure4excel["Time (min)"]-1,Sgro2015Figure4excel["Cell 1 FRET Trace (5 min pulse)"],
                               'k',linewidth=trace_width)
ax0l.plot(Sgro2015Figure4excel["Time (min)"]-1,Sgro2015Figure4excel["Cell 2 FRET Trace (5 min pulse)"],
                               'dimgrey',linewidth=trace_width)
ax0l.plot(Sgro2015Figure4excel["Time (min)"]-1,Sgro2015Figure4excel["Cell 3 FRET Trace (5 min pulse)"],
                               'darkgrey',linewidth=trace_width)
ax0l.set_xlim([0,23]); ax0l.tick_params(axis='both', which='major', labelsize=tick_font_size)
for i in range(4):
    ax0l.axvspan(4.5+i*6, 4.5+i*6+5, alpha=0.2, color=expcolor)
ax0l.text(0.5, -0.6, 'Time(min)',ha='center',va='center',
     transform = ax0l.transAxes, color = 'k', fontsize=sublabel_font_size)

#ax00= fig3.add_subplot(grid[0:2, 1])
ax00 = fig3.add_axes([0.59, 0.7, 0.31,0.18])
PeriodExp = np.linspace(3, 6, num=4); PkWdthExp = np.linspace(1,5,5)
entrainmentRs = Sgro2015Figure4["entrainmentRs"][:,0:4]; entrainmentRs[entrainmentRs == 0] = 'nan'
heatmap = ax00.pcolor(PeriodExp, PkWdthExp,np.flip(entrainmentRs,0), cmap='jet') 
heatmap.set_clim(0,1)
ax00.set_xlabel('Period(min)', size=sublabel_font_size)
ax00.set_ylabel('Peak Width(min)', size=sublabel_font_size)

cbar=fig3.colorbar(heatmap, ax=ax00,ticks=[0,0.2, 0.4,0.6,0.8,1]);
cbar.ax.tick_params(labelsize = tick_font_size) 
cbar.set_label( 'Entrainment quality',size= sublabel_font_size)
ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax00.text(-0.3 , 1.25, 'B', ha='center',va='center',
     transform = ax00.transAxes, color = 'k', fontsize=abcd_font_size)

ax1= fig3.add_subplot(grid[2:4, 0], xticks = [0.8,1.2,1.6], yticks = [0.3,0.6,0.9,1.2])
heatmap = ax1.pcolor(period_space_Gold, PkWdth_space_Gold,MeanR_Gold, cmap='jet') 
# heatmap = ax1.pcolor(period_space_Gold, PkWdth_space_Gold,MeanRnew_Gold, cmap='jet') 
heatmap.set_clim(0,1)
cbar=fig3.colorbar(heatmap, ax=ax1, ticks=[0,0.2, 0.4,0.6,0.8,1]);
# cbar=fig3.colorbar(heatmap, ax=ax1);
# cbar=fig3.colorbar(heatmap, ax=ax1, ticks=[0.96,0.98,1]);
cbar.ax.tick_params(labelsize = tick_font_size) 

ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Receptor Desensitization',color=goldcolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax1.text(-0.42 , 1.1, 'C', ha='center',va='center',
     transform = ax1.transAxes, color = 'k', fontsize=abcd_font_size)

ax2= fig3.add_subplot(grid[2:4, 1], xticks = [0.8,1.2,1.6], yticks = [0.3,0.6,0.9,1.2])
heatmap = ax2.pcolor(period_space_Maeda, PkWdth_space_Maeda,MeanR_Maeda, cmap='jet') # cmap='jet'
heatmap.set_clim(0,1)
cbar=fig3.colorbar(heatmap, ax=ax2, ticks=[0,0.2, 0.4,0.6,0.8,1]);
cbar.ax.tick_params(labelsize = tick_font_size) 

ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('CDINFB',color=maedacolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# ax2.text(-0.3 , 1.1, 'D',ha='center',va='center',
#      transform = ax2.transAxes, color = simcolor, fontsize=abcd_font_size)


ax3= fig3.add_subplot(grid[4:,0], xticks = [0.8,1.2,1.6], yticks = [0.3,0.6,0.9,1.2])
heatmap = ax3.pcolor(period_space_Sgro, PkWdth_space_Sgro,MeanR_Sgro, cmap='jet') # cmap='jet'
heatmap.set_clim(0,1)
cbar=fig3.colorbar(heatmap, ax=ax3, ticks=[0,0.2, 0.4,0.6,0.8,1]);
cbar.ax.tick_params(labelsize = tick_font_size) 

ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('IPNFB',  color= sgrocolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# ax3.text(-0.3 , 1.1, 'E', ha='center',va='center',
#      transform = ax3.transAxes, color = simcolor, fontsize=abcd_font_size)

ax4= fig3.add_subplot(grid[4:,1], xticks = [0.8,1.2,1.6], yticks = [0.3,0.6,0.9,1.2])
heatmap = ax4.pcolor(period_space_Kamino, PkWdth_space_Kamino,MeanR_Kamino, cmap='jet') # cmap='jet'
heatmap.set_clim(0,1)
cbar=fig3.colorbar(heatmap, ax=ax4, ticks=[0,0.2, 0.4,0.6,0.8,1]);
cbar.ax.tick_params(labelsize = tick_font_size) 
#ax4.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
#ax4.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('IFFL', color= kaminocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# ax4.text(-0.3 , 1.1, 'F', ha='center',va='center',
#      transform = ax4.transAxes, color = simcolor, fontsize=abcd_font_size)

fig3.text(0.5, 0.05,'Period, A.U.', ha='center', va='center',fontsize=label_font_size)
fig3.text(0.04, 0.4, 'Peak Width, A.U.', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

# fig3.savefig('fig4_sc_entrainment_200530.png', bbox_inches='tight')
# fig3.savefig('fig4_sc_entrainment_goldbeter0_1_tight_200530.png', bbox_inches='tight')
#%% fig S5 Gldbeter entrainment
from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_agent_3var

abcd_font_size = 28
label_font_size=22
title_font_size = 24
sublabel_font_size = 18
trace_width=3
tick_font_size=20

k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10; L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01; theta=0.01
e= 1
q=4000
sig=0.6 # compared to 0.57
v=12; k= 4 # k prime in the paper
ki=1.7 
kt=0.9
kc=5.4 
h=5
Goldbeter3AgentParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}
p0=0.8; a0=3; b0=0.9; g0=0
dt=0.001 ; # t_tot=20*Nt; t=list(np.arange(0,t_tot,dt))
cAMP = 1 # extracellular cAMP
NumofCycle = 8 # cycle stimulation

period_space = [1];PkWdth_space = [0.3,0.6,0.9]
count=1

fig3 = plt.figure(figsize=(cm2inch(38),cm2inch(13))); 
grid = plt.GridSpec(3, 9, wspace=0.1, hspace=1.6)

ax0= fig3.add_subplot(grid[0:2, 0:2], xticks = [0.8,1.2,1.6], yticks = [0.3,0.6,0.9,1.2])
heatmap = ax0.pcolor(period_space_Gold, PkWdth_space_Gold,MeanR_Gold, cmap='jet') 
cbar=fig3.colorbar(heatmap, ax=ax0, ticks=[0.96,0.98,1]);
cbar.ax.tick_params(labelsize = tick_font_size) 
cbar.set_label( 'Entrainment quality',size= sublabel_font_size)

ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_title('Receptor desensitization',color=goldcolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax0.text(-0.42 , 1.16, 'A', ha='center',va='center',
     transform = ax0.transAxes, color ='k', fontsize=abcd_font_size)

fig3.text(0.34, 0.95,'B', ha='center', va='center',fontsize = abcd_font_size) 

fig3.text(0.64, 0,'Time, A.U.', ha='center', va='center',fontsize = sublabel_font_size) 
fig3.text(0.35, 0.81,r'$cAMP_{e}$'+'\ninput',  ha='center', va='center',
          rotation='vertical', fontsize = sublabel_font_size) 
fig3.text(0.35, 0.32,r'$cAMP_{i}$ response',  ha='center', va='center',
          rotation='vertical',fontsize = sublabel_font_size) 

for period in period_space: #j_test:# range(len(period_space_Gold)):# 
    # period = period_space_Gold[j]
    for PkWdth in PkWdth_space: #k in k_test:#range(len(PkWdth_space_Gold)): #
        # PkWdth = PkWdth_space_Gold[k]
        if PkWdth>=period: 
            pass
        else:
            stim_1cycle=np.zeros(int(period*Nt_Goldbeter/dt))
            stim_1cycle[0:int(PkWdth*Nt_Goldbeter/dt)]=cAMP;
            signal_trace=np.concatenate((np.zeros(int(1*Nt_Goldbeter/dt)), np.tile(stim_1cycle,(NumofCycle)),np.zeros(int(1*Nt_Goldbeter/dt))),axis=0);
            
            t = list(np.arange(0,len(signal_trace)*dt,dt))
            t_plot_Gold=np.array(t)/Nt_Goldbeter
            
            p_trace=[p0]; b_trace=[b0]; g_trace=[g0]
            Goldbeter3_agent=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
            for i in range(len(t)-1):
                p_now=p_trace[i]
                p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,signal_trace[i])
                p_trace.append(p_next)
                b_trace.append(b_next)
                g_trace.append(g_next)
                      
            b_trace = np.array(b_trace); b_trace = b_trace/Nh_Goldbeter

            r=np.zeros(NumofCycle-1) # list that stores correlation coefficient to the first peak
            r_new = np.zeros(NumofCycle-1) # correlation coeff take in height difference
            InitPk=b_trace[int(1*Nt_Goldbeter/dt): int((1+period)*Nt_Goldbeter/dt)]
            for m in range(NumofCycle-1):
                FollowPk=b_trace[int((1+(m+1)*period)*Nt_Goldbeter/dt) : (int((1+(m+1)*period)*Nt_Goldbeter/dt)+len(InitPk))]
                R = np.corrcoef(InitPk,FollowPk)
                r[m]=R[1,0]
                r_new[m]=r[m]*(1-(math.sqrt(np.var(InitPk))-math.sqrt(np.var(FollowPk)))/math.sqrt(np.var(InitPk)))
                      
            # plot supplimentary figure         
            ax1= fig3.add_subplot(grid[0, count*2+1:count*2+3])
            ax1.plot(t_plot_Gold,signal_trace,linewidth = trace_width-1, color = goldcolor)
            # if count ==0:
            #     ax1.set_ylabel(r'$cAMP_{e}$ input', fontsize= label_font_size); 
            ax1.set_title('Period '+str(period)+', stimulation\npeak width '+ str(PkWdth), fontsize = sublabel_font_size)
            ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
            ax1.set_xlim([0, 10])
            ax2= fig3.add_subplot(grid[1:, count*2+1:count*2+3])
            ax2.plot(t_plot_Gold,b_trace,linewidth = trace_width, color = goldcolor)
            # if count ==0:
            #     ax2.set_ylabel('Entrainment quality', fontsize= label_font_size);
            ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
            ax2.set_title('Entrainment\n quality= {:.4f}'.format(np.mean(r)), fontsize = sublabel_font_size)
            ax2.set_xlim([0, 10])
            count=count+1
            plt.show()
grid.tight_layout(fig3,rect=[0, 0, 1, 1])

fig3.savefig('figS4_Goldbeter_entrainment_tight_200530.png', bbox_inches='tight')

#%% Fig 5 single cell FCD
# load simulation data
import_npz('../model_outputs/single_cell_FCD_200322.npz',globals())
# Load experimental data
my_dir = r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/'
Kamino_FCD = pd.read_excel(my_dir+r'Kamino_FCD_exp_data.xlsx',sheet_name='Sheet1')

# Panel B: demonstration of concepts
from Kamino2017_agent_and_pop_FUN import Kamino2017_agent 
tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
AgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
x0=0.01; y0=0.05; z0=0.005
Kamino_agent=Kamino2017_agent([x0,y0,z0],AgentParam)
dt=0.001; t_tot=30 * Nt_Kamino; t=list(np.arange(0,t_tot,dt))

z0First = 1; FC = 3 
signal_trace=np.zeros(len(t))
stim_time_step1=int(round(0.33*t_tot/dt)) ; signal_trace[stim_time_step1:] = z0First
stim_time_step2=int(round(0.66*t_tot/dt)) ; signal_trace[stim_time_step2:] = FC*z0First
x_trace_B=[x0]; y_trace_B=[y0]
        
for i in range(len(t)-1):
    x_now=x_trace_B[i]
    y_now=y_trace_B[i]
    x_next,y_next,z_next= Kamino_agent.update( dt, signal_trace[i])
    x_trace_B.append(x_next)
    y_trace_B.append(y_next)
# Convert into np array
x_trace_B = np.array(x_trace_B) # vectorize p_trace
y_trace_B = np.array(y_trace_B)
t_plot_Kamino_B = np.array(t)/Nt_Kamino

# set up plotting specs
abcd_font_size = 28
title_font_size = 19
label_font_size = 17 # 16
sublabel_font_size = 17
tick_font_size = 17
legend_font_size = 14
trace_width = 3
annotate_size = 15
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Gold)+1))

fig3 = plt.figure(figsize=(cm2inch(21),cm2inch(26))) 
grid = plt.GridSpec(7, 2, wspace=0.25, hspace=1.25)

ax0u = fig3.add_axes([0.12, 0.88, 0.315,0.075],xticks=[0,10,20,30],yticks=[0,2])
ax0u.plot(t_plot_Kamino_B,signal_trace, 'k',linewidth=trace_width)
ax0u.set_xlim([0,30]); ax0u.set_ylim([-0.5,3.5]); 
ax0u.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0u.text(-0.33 , 1.3, 'A', ha='center',va='center',
     transform = ax0u.transAxes, color = expcolor, fontsize=abcd_font_size)

ax0u.text(-0.28, 0.5, r'$cAMP_{e}$'+'\n input',ha='center',va='center', rotation='vertical',
     transform = ax0u.transAxes, color = 'k', fontsize=sublabel_font_size)

ax0u.annotate(s='', xy=(11,-0.5), xytext=(11,1.2), arrowprops={'arrowstyle':'<->', 'lw':2.5})
ax0u.text(0.5, 0.65, 'Prime\nconc.',ha='center',va='center', 
     transform = ax0u.transAxes, color = 'k', fontsize=annotate_size)

ax0u.annotate(s='', xy=(21,0.7), xytext=(21,3.2), arrowprops={'arrowstyle':'<->', 'lw':2.5})
ax0u.text(0.86, 0.46, 'Fold\nchange',ha='center',va='center', 
     transform = ax0u.transAxes, color = 'k', fontsize=annotate_size)

ax0l = fig3.add_axes([0.12, 0.77, 0.315,0.075],xticks=[0,10,20,30],yticks=[0,0.2,0.4])
ax0l.plot(t_plot_Kamino_B, y_trace_B, linewidth=trace_width)
ax0l.set_xlim([0,30]); ax0l.set_ylim([-0.05,0.4]); 
ax0l.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0l.set_xlabel('Time, A.U.', size=sublabel_font_size)
ax0l.text(-0.28, 0.5, r'$cAMP_{i}$'+'\n response',ha='center',va='center', rotation='vertical',
     transform = ax0l.transAxes, color = 'k', fontsize=sublabel_font_size)

ax0l.annotate('', xy=(23,0.01), xytext=(23,0.2),rotation = 'vertical',
            arrowprops={'arrowstyle': '<->','lw':2.5}, va='center')
ax0l.text(0.75, 0.72, 'Second peak\nprominence',ha='center',va='center', 
      transform = ax0l.transAxes, color = 'k', fontsize=annotate_size)

ax0 = fig3.add_axes([0.6, 0.77, 0.32,0.18],yticks=[1.5,2,2.5])
#for i in range(len(z0First_space_Sgro)):
#    ax3.plot(FC_space_Sgro_noise,PkPrm_Sgro_mean_noise[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Sgro_noise[i]))
#plt.errorbar(Kamino_FCD["FC_100pM"], Kamino_FCD["100pM mean"], yerr=None, xerr=None, fmt='', ecolor=None, elinewidth=None, capsize=None, barsabove=False, lolims=False, uplims=False, xlolims=False, xuplims=False, errorevery=1, capthick=None, *, data=None, **kwargs)
ax0.errorbar(Kamino_FCD["FC_100pM"], Kamino_FCD["100pM mean"], yerr=Kamino_FCD["100pM SD"], xerr=None, color='SteelBlue', linewidth=trace_width, label='100pM', ecolor='SteelBlue', elinewidth=trace_width,capsize=5,capthick=2)
ax0.errorbar(Kamino_FCD["FC_1nM"], Kamino_FCD["1nM mean"], yerr=Kamino_FCD["1nM SD"], xerr=None, color='SkyBlue', linewidth=trace_width, label='1nM', ecolor='SkyBlue', elinewidth=trace_width,capsize=5,capthick=2)
ax0.errorbar(Kamino_FCD["FC_3nM"], Kamino_FCD["3nM mean"], yerr=Kamino_FCD["3nM SD"], xerr=None, color='turquoise', linewidth=trace_width, label='3nM', ecolor='turquoise', elinewidth=trace_width,capsize=5,capthick=2)
ax0.errorbar(Kamino_FCD["FC_10nM"], Kamino_FCD["10nM mean"], yerr=Kamino_FCD["10nM SD"], xerr=None, color='cyan', linewidth=trace_width, label='10nM', ecolor='cyan', elinewidth=trace_width,capsize=5,capthick=2)
ax0.set_ylim([1.5,2.8])
# ax0.set_ylabel( 'Response Amplitude, A.U.',fontsize=tick_font_size)
ax0.set_xlabel(r'$cAMP_{e}$'+' fold change',fontsize=label_font_size)
ax0.set_xscale('log')
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_title('Experiment', color = 'k',fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax0.legend();
ax0.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})
#ax0.text(-0.1, 1.15, 'B', ha='center',va='center',
#         transform = ax0.transAxes, color = 'b', fontsize=abcd_font_size)
ax0.text(-0.27, 0.53, 'Second peak \n prominence, A.U.', ha='center',va='center', rotation ='vertical',
         transform = ax0.transAxes, fontsize=label_font_size)
#
ax1= fig3.add_subplot(grid[2:4, 0])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Gold)+1) )
for i in range(len(z0First_space_Gold)):
    ax1.plot(FC_space_Gold,PkPrm_Gold[i,:], color=colors[i],
             linewidth=trace_width,label='Priming Conc. '+str(int(z0First_space_Gold[i]*10))+' , A.U.')
# ax1.set_ylim([-0.1,1.2]); 
ax1.set_xlim([FC_space_Gold[0],FC_space_Gold[-1]])
#ax1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#ax1.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax1.set_xscale('log')
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
# ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax1.set_title('Receptor desensitization', color= goldcolor,
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax1.text(-0.3, 1.15, 'B', ha='center',va='center',
         transform = ax1.transAxes, color = simcolor, fontsize=abcd_font_size)
#leg = ax1.legend();
#ax1.legend( frameon=False,loc='bottom center',ncol=1,prop={'size': legend_font_size})

ax2= fig3.add_subplot(grid[2:4, 1])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Loomis)+1) )
for i in range(len(z0First_space_Loomis)):
    ax2.plot(FC_space_Loomis,PkPrm_Loomis[i,:],color=colors[i], linewidth=trace_width,label='Prime Conc.'+str(z0First_space_Loomis[i]))

# ax2.set_ylim([0,0.6]); 
ax2.set_xlim([FC_space_Loomis[0],FC_space_Loomis[-1]])
ax2.set_xscale('log')
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('CDINFB', color = maedacolor,
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# ax2.text(-0.1, 1.18, 'C', ha='center',va='center',
#          transform = ax2.transAxes, color = simcolor, fontsize=abcd_font_size)

# Sgro with noise
ax5= fig3.add_subplot(grid[4:6, 0])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Sgro)+1) )
for i in range(len(z0First_space_Sgro_noise)):
    ax5.plot(FC_space_Sgro_noise, PkPrm_Sgro_mean_noise[i,:],'o-', color = colors[i], lw = trace_width, 
             ms = 4,  label='Prime Conc.'+str(z0First_space_Sgro_noise[i]))
    ax5.errorbar(FC_space_Sgro_noise, PkPrm_Sgro_mean_noise[i,:], yerr=PkPrm_Sgro_se_noise[i,:],
                 fmt = 'o', color=colors[i], ecolor= colors[i], elinewidth=trace_width, capsize=5, capthick=2)
ax5.set_ylim([-0.1,1]);
ax5.set_xlim([FC_space_Sgro[0],FC_space_Sgro[-1]])
ax5.set_xscale('log')
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5.set_title('IPNFB', color=sgrocolor,
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# ax5.text(-0.1, 1.18, 'D', ha='center',va='center',
#          transform = ax5.transAxes, color = simcolor, fontsize=abcd_font_size)

ax4= fig3.add_subplot(grid[4:6, 1])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Kamino)+1) )
for i in range(len(z0First_space_Kamino)):
    ax4.plot(FC_space_Kamino,PkPrm_Kamino[i,:],color=colors[i], linewidth=trace_width,label='Prime Conc.'+str(z0First_space_Kamino[i]))
ax4.set_ylim([0,1.1]); 
ax4.set_xlim([FC_space_Kamino[0],FC_space_Kamino[-1]])    
ax4.set_xscale('log')
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('IFFL',  color=kaminocolor,
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# ax4.text(-0.1, 1.18, 'E', ha='center',va='center',
#          transform = ax4.transAxes, color = simcolor, fontsize=abcd_font_size)
#leg = ax4.legend();
#ax4.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})
ax4.legend(bbox_to_anchor=(0.6,-0.32),ncol=2,prop={'size': 16})

fig3.text(0.5, 0.18, r'$cAMP_{e}$'+' fold change', ha='center', va='center',fontsize=label_font_size)
fig3.text(0.05, 0.45, 'Second peak prominence', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

# plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
fig3.savefig('fig5_sc_FCD_tight_200530.png', bbox_inches='tight') 

#%% Fig S5 FCD second input concentration
# load simulation data
import_npz('../model_outputs/figS5_sc_FCD_2ndconc_042320.npz',globals())

# set up plotting specs
abcd_font_size = 28
label_font_size=22
title_font_size = 24
sublabel_font_size = 18
trace_width=3
tick_font_size=20

colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Gold)+1))

fig3 = plt.figure(figsize=(cm2inch(22),cm2inch(15))) 
grid = plt.GridSpec(5, 2, wspace=0.2, hspace=0.5)

ax1= fig3.add_subplot(grid[0:4, 0])
for i in range(len(z0First_space_Gold)):
    
    ax1.plot(scd_space_Gold,PkPrm_Gold[i,:], linewidth=trace_width, color = colors[i],
             label='Primie conc.='+str(z0First_space_Gold[i])+ ' A.U.')

ax1.set_ylabel( 'second spike prominence',fontsize=label_font_size)
ax1.set_xscale('log')
#ax1.set_xlabel('Fold change in extracellular cAMP',fontsize=label_font_size)
ax1.set_xlim([scd_space_Gold[0],scd_space_Gold[-1]])
# ax1.set_xlabel(r'Second step [$cAMP_{e}$]',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Receptor desensitization', color = goldcolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax2= fig3.add_subplot(grid[0:4, 1])
for i in range(len(z0First_space_Loomis)):
    ax2.plot(scd_space_Loomis,PkPrm_Loomis[i,:], linewidth=trace_width, color = colors[i],
             label='Prime conc.= '+str(z0First_space_Loomis[i]))
ax2.set_xscale('log')
ax2.set_xlim([scd_space_Loomis[0],scd_space_Loomis[-1]])
# ax2.set_xlabel(r'Fold change in [$cAMP_{e}$]',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('CDINFB', color = maedacolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax2.legend(bbox_to_anchor=(0.74,-0.2),ncol=2,prop={'size': 18})

fig3.text(0.5, 0.18,'Second step '+ r'$cAMP_{e}$', 
          ha='center', va='center',fontsize=label_font_size)

ax1.text(-0.2, 1.1, 'A', ha='center',va='center',
         transform = ax1.transAxes, color = 'k', fontsize=abcd_font_size)
ax2.text(-0.2, 1.1, 'B', ha='center',va='center',
         transform = ax2.transAxes, color = 'k', fontsize=abcd_font_size)
plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)

fig3.savefig('figS5_sc_FCD_2ndinput_tight_200530.png', bbox_inches='tight') 

#%% Figure S6 FCD Sgro violin plot

colors = plt.cm.Greens(np.linspace(0,1,len(z0First_space_Sgro_noise)+2))

title_font_size = 20
label_font_size = 18
tick_font_size = 17
legend_font_size = 17
trace_width = 5
markers=['+', 'd', '2', 'x']
# plot FC vs. second response amplitude, violoin plot
fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(cm2inch(45),cm2inch(13)))
for i in range(len(z0First_space_Sgro_noise)):
    z0_now = z0First_space_Sgro_noise[i]
    data = PkPrm_Sgro_noise[i,:,:].tolist()
    violin_parts = ax[i].violinplot(data,showmeans=True, showmedians=False, showextrema=False)
    v_mean = violin_parts['cmeans']
    v_mean.set_edgecolor('k')
    v_mean.set_linewidth(3)
    for vp in violin_parts['bodies']:
        vp.set_facecolor(colors[i+2])
        vp.set_edgecolor(colors[i+2])
        vp.set_alpha(0.8)    
    ax[i].set_xticks([1,2,3,4,5,6,7,8]); 
    ax[i].set_xticklabels(np.around(FC_space_Sgro_noise,decimals=1), 
                          rotation=45,fontsize=tick_font_size)
    ax[i].tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax[i].set_title('Prime concentration '+str(round(z0_now,2)), fontdict={'fontsize': label_font_size, 'fontweight': 'medium'})
fig.text(0.5, 0.97, 'IPNFB (with noise)', ha='center', va='center',fontsize=title_font_size,color=sgrocolor,)
fig.text(0.5, 0.02, r'$cAMP_{e}$'+' fold change', ha='center', va='center',fontsize=label_font_size)
fig.text(0.01, 0.55, 'Second peak prominence, A.U.', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

plt.subplots_adjust(left=0.05, right=0.95, top=0.85, bottom=0.18, wspace=0.15, hspace=0)
fig.savefig('figS6_sc_FCD_Sgro_violin_tight_200530.png', bbox_inches='tight') 

#%% Fig 6 pop oscillations
# experimental data
my_dir = r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/'
Sgro2015Figure6excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure6')

# import simulation data
import_npz('../model_outputs/pop_oscillation_SCnoise_200409.npz',globals())

abcd_font_size = 28
label_font_size=22
title_font_size = 22
sublabel_font_size = 22
trace_width=3
tick_font_size=20
SC_traces_idx = [0,2,4,6,8,10]

fig3 = plt.figure(figsize=(cm2inch(27),cm2inch(19)))
grid = plt.GridSpec(3, 2, wspace=0.15, hspace=0.78)

ax0= fig3.add_subplot(grid[0, 0])
ax0.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["No External cAMP Mean Trace"],
                              color = 'k',linewidth=trace_width)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_title('Experiment', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax0.set_ylim([-0.1,0.5]); ax0.set_xlim([0,120])
ax0.set_xlabel('Time (min)', size=sublabel_font_size)
ax0.text(-0.14 , 1.2, 'A', ha='center',va='center',
     transform = ax0.transAxes, color = 'k', fontsize=abcd_font_size)

ax1= fig3.add_subplot(grid[0, 1])
for this_idx in SC_traces_idx:
    this_trace = b_trace_norm_later[this_idx,:] 
    ax1.plot(t_plot_Goldbeter_short,this_trace, color='k',alpha=0.6, linewidth=2)     
# Plot population mean
ax1.plot(t_plot_Goldbeter_short,b_trace_norm_mean_later,
         color=goldcolor,alpha=0.8, linewidth=trace_width)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Receptor Desensitization', color=goldcolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax1.set_ylim([-0.3,1.7]); ax1.set_xlim([0,15])
ax1.text(-0.08 , 1.2, 'B', ha='center',va='center',
     transform = ax1.transAxes, color = 'k', fontsize=abcd_font_size)

ax2= fig3.add_subplot(grid[1,0])
for this_idx in SC_traces_idx:
    this_trace = cAMPi_trace_norm_later[this_idx,:] 
    ax2.plot(t_plot_Maeda_short,this_trace, color='k',alpha=0.6, linewidth=2)     
# Plot population mean
ax2.plot(t_plot_Maeda_short,cAMPi_trace_norm_mean_later,
         color=maedacolor,alpha=0.8, linewidth=trace_width)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('CDINFB',color= maedacolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax2.set_ylim([0.1,0.6]); ax2.set_xlim([0,15])
ax2.text(-0.14 , 1.2, 'C', ha='center',va='center',
     transform = ax2.transAxes, color = 'k', fontsize=abcd_font_size)

ax3= fig3.add_subplot(grid[1, 1])
for this_idx in SC_traces_idx:
    this_trace = gregor_campCyto_trace_norm_later[this_idx,:] 
    ax3.plot(t_plot_Gregor_short,this_trace, color='k',alpha=0.6, linewidth=2)     
# Plot population mean
ax3.plot(t_plot_Gregor_short,gregor_campCyto_trace_norm_mean_later,
         color=gregorcolor,alpha=0.8, linewidth=trace_width)
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Phase Oscillator', color=gregorcolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax3.set_ylim([-0.2,1.5]); ax3.set_xlim([0,15])
ax3.text(-0.08 , 1.2, 'D', ha='center',va='center',
     transform = ax3.transAxes, color = 'k', fontsize=abcd_font_size)

ax4= fig3.add_subplot(grid[2, 0])
for this_idx in SC_traces_idx:
    this_trace = A_trace_norm_later[this_idx,:] 
    ax4.plot(t_plot_Sgro_short,this_trace, color='k',alpha=0.6, linewidth=2)     
# Plot population mean
ax4.plot(t_plot_Sgro_short,A_trace_norm_mean_later,
         color= sgrocolor,alpha=0.8, linewidth=trace_width)
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('IPNFB', color=sgrocolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax4.set_ylim([-0.3,1.3]); ax4.set_xlim([0,15])
ax4.text(-0.14 , 1.2, 'E', ha='center',va='center',
     transform = ax4.transAxes, color = 'k', fontsize=abcd_font_size)

ax5= fig3.add_subplot(grid[2,1])
for this_idx in SC_traces_idx:
    this_trace = y_trace_norm_later[this_idx,:] 
    ax5.plot(t_plot_Kamino_short,this_trace, color='k',alpha=0.6, linewidth=2)     
# Plot population mean
ax5.plot(t_plot_Kamino_short,y_trace_norm_mean_later,
         color= kaminocolor,alpha=0.8, linewidth=trace_width)
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5.set_title('IFFL', color= kaminocolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax5.set_ylim([-0.5,1.5]); ax5.set_xlim([0,15])
ax5.text(-0.08 , 1.2, 'F', ha='center',va='center',
     transform = ax5.transAxes, color = 'k', fontsize=abcd_font_size)

fig3.text(0.04, 0.5, r'$cAMP_{i}$', fontsize=label_font_size,va='center', rotation='vertical')
fig3.text(0.56, 0.02, 'Time, A.U.', fontsize=label_font_size, ha='center')
plt.show()

plt.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.1)
# grid.tight_layout(fig3,rect=[0, 0, 1, 1])
# fig3.savefig('fig6_pop_oscillation_tight_200530.png', bbox_inches='tight') 

#%% Fig 7 population add cAMP, with single cell noise
# Experimental data
my_dir = r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/'
Sgro2015Figure6excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',
                                     sheet_name='Figure6')
# load saved npz output file 
import_npz('../model_outputs/Fig8_pop_add_cAMP_042320.npz',globals())
#%%
title_font_size = 20
label_font_size = 20
sublabel_font_size = 18
tick_font_size = 16
text_size = 14
trace_width = 2
abcd_font_size = 28

fig = plt.figure(figsize=(cm2inch(45),cm2inch(34)))
grid = plt.GridSpec(6, 3, wspace=0.17, hspace=0.8)
# experiment
ax01= fig.add_subplot(grid[0, 0])
ax01.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["Low External cAMP Mean Trace"],
                               color = 'k', linewidth=trace_width)
ax01.axvspan(60, 120, alpha=0.2, color='b')
ax01.set_ylim([-0.1,0.6]); ax01.set_xlim([0,120])
ax01.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax01.text(0.7,0.75,' Low External cAMP, \n 5-10nM', horizontalalignment='center',verticalalignment='center',
         transform = ax01.transAxes, color = 'k', fontsize= text_size)
ax01.set_title('Experiment',size=title_font_size)
ax02= fig.add_subplot(grid[1, 0])
ax02.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["Intermediate External cAMP Mean Trace"],
                               color = 'k', linewidth=trace_width)
ax02.axvspan(60, 120, alpha=0.2, color = expcolor)
ax02.set_ylim([-0.1,0.6]); ax02.set_xlim([0,120])
ax02.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax02.text(0.7,0.75,' Intermediate External \n cAMP, 10-20nM', horizontalalignment='center',verticalalignment='center',
         transform = ax02.transAxes, color = 'k', fontsize=text_size)
ax02.set_ylabel('FRET Signal, A.U.', size=sublabel_font_size)

ax03= fig.add_subplot(grid[2, 0])
ax03.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["High External cAMP Mean Trace"],
                               color = 'k', linewidth=trace_width)
ax03.axvspan(60, 120, alpha=0.2, color='b')
ax03.set_ylim([-0.1,0.6]);  ax03.set_xlim([0,120])
ax03.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax03.text(0.7,0.75,' High External cAMP, \n 100nM', horizontalalignment='center',verticalalignment='center',
         transform = ax03.transAxes, color = 'k', fontsize=text_size)
# ax03.set_xlabel('Time (min)', size=sublabel_font_size)
ax01.text(-0.06 , 1.2, 'A',
         horizontalalignment='center',verticalalignment='center',
         transform = ax01.transAxes, color = expcolor, fontsize=abcd_font_size)
# Goldbeter 1987
ax11= fig.add_subplot(grid[0, 1])
for count in range(5):
    ax11.plot(t_plot_Goldbeter,b_traces_norm[0,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax11.plot(t_plot_Goldbeter,b_traces_norm_norm_mean[0,:], color= goldcolor,linewidth=trace_width)
ax11.text(0.73,0.88,r'Low $cAMP_{e}$ input', ha='center',va='center',
     transform = ax11.transAxes, color = 'k', fontsize= text_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax11.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax11.axvspan(15, 30, alpha=0.2, color= simcolor)
ax11.set_xlim([0,30]); ax11.set_ylim([-0.25,1.75])
ax11.text(-0.06 , 1.2, 'B',
         horizontalalignment='center',verticalalignment='center',
         transform = ax11.transAxes, color = simcolor, fontsize=abcd_font_size)
ax11.set_title('Receptor desensitization',color = goldcolor, size=title_font_size)

ax12= fig.add_subplot(grid[1, 1])
for count in range(5):
    ax12.plot(t_plot_Goldbeter, b_traces_norm[1,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax12.plot(t_plot_Goldbeter, b_traces_norm_norm_mean[1,:], color=goldcolor,linewidth=trace_width)
    

ax12.text(0.73,0.88,r'Intermediate $cAMP_{e}$'+' input', ha='center',va='center',
     transform = ax12.transAxes, color = 'k', fontsize= text_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax12.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax12.axvspan(15, 30, alpha=0.2, color=simcolor)
ax12.set_xlim([0,30]); ax12.set_ylim([-0.25,1.75])

ax13= fig.add_subplot(grid[2, 1])
for count in range(5):
    ax13.plot(t_plot_Goldbeter, b_traces_norm[2,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax13.plot(t_plot_Goldbeter, b_traces_norm_norm_mean[2,:], color=goldcolor,linewidth=trace_width)
ax13.text(0.73,0.88,r'High $cAMP_{e}$ input', ha='center',va='center',
     transform = ax13.transAxes, color = 'k', fontsize= text_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax13.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax13.axvspan(15, 30, alpha=0.2, color=simcolor)
ax13.set_xlim([0,30]); ax13.set_ylim([-0.25,1.75])
# CDINFB
ax21= fig.add_subplot(grid[0, 2])
for count in range(5):
    ax21.plot(t_plot_Maeda_short,cAMPi_traces_norm[0,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax21.plot(t_plot_Maeda_short, cAMPi_traces_norm_mean[0,:], color=maedacolor,linewidth=trace_width)
ax21.text(0.73,0.88,r'Low $cAMP_{e}$ input', ha='center',va='center',
     transform = ax21.transAxes, color = 'k', fontsize= text_size)
ax21.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax21.axvspan(15, 30, alpha=0.2, color = simcolor)
ax21.set_xlim([0,30]); ax21.set_ylim([0.1,0.9])
ax21.text(-0.06 , 1.2, 'C',
         horizontalalignment='center',verticalalignment='center',
         transform = ax21.transAxes, color = simcolor, fontsize=abcd_font_size)
ax21.set_title('CDINFB',color = maedacolor, size=title_font_size)

ax22= fig.add_subplot(grid[1, 2])
for count in range(5):
    ax22.plot(t_plot_Maeda_short,cAMPi_traces_norm[1,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax22.plot(t_plot_Maeda_short, cAMPi_traces_norm_mean[1,:], 
          color= maedacolor,linewidth=trace_width)
ax22.text(0.73,0.88,r'Intermediate $cAMP_{e}$ input', ha='center',va='center',
          transform = ax22.transAxes, color = 'k', fontsize= text_size)
ax22.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax22.axvspan(15, 30, alpha=0.2, color=simcolor)
ax22.set_xlim([0,30]); ax22.set_ylim([0.1,0.9])

ax23= fig.add_subplot(grid[2, 2])
for count in range(5):
    ax23.plot(t_plot_Maeda_short,cAMPi_traces_norm[2,count,:],
             color= 'darkgrey',alpha=0.5, linewidth=trace_width-1)
ax23.plot(t_plot_Maeda_short, cAMPi_traces_norm_mean[2,:], 
          color=maedacolor,linewidth=trace_width)
ax23.text(0.73,0.88,r'High $cAMP_{e}$ input', ha='center',va='center',
     transform = ax23.transAxes, color = 'k', fontsize= text_size)
ax23.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax23.axvspan(15, 30, alpha=0.2, color = simcolor)
ax23.set_xlim([0,30]); ax23.set_ylim([0.1,0.9])

# Phase oscillator
ax31= fig.add_subplot(grid[3, 0])
ax31.plot(t_plot_Gregor[1:],campCyto_traces[0,1:] , color=gregorcolor,linewidth=trace_width)  
for count in range(5):
    # campCyto_traces_single_cell[0,count,:] = campCyto_traces_single_cell[0,count,:]/np.amax(campCyto_traces_single_cell[0,count,:])
    ax31.plot(t_plot_Gregor,campCyto_traces_single_cell[0,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
 
ax31.text(0.73,0.88,r'Low $cAMP_{e}$ input', ha='center',va='center',
     transform = ax31.transAxes, color = 'k', fontsize= text_size)
ax31.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax31.axvspan(15, 30, alpha=0.2, color= simcolor)
ax31.set_xlim([0,30]); ax31.set_ylim([-0.2,1.5])
ax31.text(-0.06 , 1.2, 'D', ha='center',va='center',
         transform = ax31.transAxes, color = simcolor, fontsize=abcd_font_size)
ax31.set_title('Phase oscillator',color = gregorcolor, size=title_font_size)

ax32= fig.add_subplot(grid[4, 0])
ax32.plot(t_plot_Gregor[1:],campCyto_traces[1,1:], color=gregorcolor,linewidth=trace_width)
for count in range(5):
    # campCyto_traces_single_cell[1,count,:] = campCyto_traces_single_cell[1,count,:]/np.amax(campCyto_traces_single_cell[1,count,:])
    ax32.plot(t_plot_Gregor,campCyto_traces_single_cell[1,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax32.text(0.73,0.88,r'Intermediate $cAMP_{e}$'+' input', ha='center',va='center',
     transform = ax32.transAxes, color = 'k', fontsize= text_size)
ax32.tick_params(grid_linewidth =tick_font_size, labelsize = tick_font_size)
ax32.axvspan(15, 30, alpha=0.2, color = simcolor)
ax32.set_xlim([0,30]); ax32.set_ylim([-0.2, 1.5])

ax33= fig.add_subplot(grid[5, 0])
ax33.plot(t_plot_Gregor[1:],campCyto_traces[2,1:], color= gregorcolor, linewidth=trace_width)
for count in range(5):
    # campCyto_traces_single_cell[2,count,:] = campCyto_traces_single_cell[2,count,:]/np.amax(campCyto_traces_single_cell[2,count,:])
    ax33.plot(t_plot_Gregor,campCyto_traces_single_cell[2,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax33.text(0.73,0.88,r'High $cAMP_{e}$ input', ha='center',va='center',
     transform = ax33.transAxes, color = 'k', fontsize= text_size)
ax33.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax33.axvspan(15, 30, alpha=0.2, color=simcolor)
ax33.set_xlim([0,30]); ax33.set_ylim([-0.2,1.5])

# IPNFB
ax41= fig.add_subplot(grid[3, 1])
for count in range(5):
    ax41.plot(t_plot_Sgro,A_traces_single_cell[0,count,:],
             color= 'darkgrey' ,alpha=0.3, linewidth=trace_width-1)
ax41.plot(t_plot_Sgro,A_traces[0,:], color= sgrocolor,linewidth=trace_width)
ax41.text(0.73,0.88,r'Low $cAMP_{e}$ input', ha='center',va='center',
     transform = ax41.transAxes, color = 'k', fontsize= text_size)
ax41.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax41.axvspan(15, 30, alpha=0.2, color=simcolor)
ax41.set_xlim([0,30]); ax41.set_ylim([-0.25,1.25])
ax41.text(-0.06 , 1.2, 'E',ha='center',va='center',
         transform = ax41.transAxes, color = simcolor, fontsize=abcd_font_size)
ax41.set_title('IPNFB',color = sgrocolor, size=title_font_size)

ax42= fig.add_subplot(grid[4, 1])
for count in range(5):
    ax42.plot(t_plot_Sgro,A_traces_single_cell[1,count,:],
             color='darkgrey',alpha=0.3, linewidth=trace_width-1)
ax42.plot(t_plot_Sgro,A_traces[1,:], color=sgrocolor,linewidth=trace_width)
    

ax42.text(0.73,0.88,r'Intermediate $cAMP_{e}$'+' input', ha='center',va='center',
     transform = ax42.transAxes, color = 'k', fontsize= text_size)
ax42.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax42.axvspan(15, 30, alpha=0.2, color=simcolor)
ax42.set_xlim([0,30]); ax42.set_ylim([-0.25,1.25])

ax43= fig.add_subplot(grid[5, 1])
for count in range(5):
    ax43.plot(t_plot_Sgro,A_traces_single_cell[2,count,:],
             color='darkgrey',alpha=0.3, linewidth=trace_width-1)
ax43.plot(t_plot_Sgro,A_traces[2,:], color=sgrocolor,linewidth=trace_width)
ax43.text(0.73,0.88,r'High $cAMP_{e}$ input', ha='center',va='center',
     transform = ax43.transAxes, color = 'k', fontsize= text_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax43.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax43.axvspan(15, 30, alpha=0.2, color=simcolor)
ax43.set_xlim([0,30]); ax43.set_ylim([-0.25,1.25])

# IFFL
ax51= fig.add_subplot(grid[3, 2])
for count in range(5):
    ax51.plot(t_plot_Kamino,y_traces_norm[0,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax51.plot(t_plot_Kamino,y_traces_norm_mean[0,:], color= kaminocolor,linewidth=trace_width)
ax51.text(0.73,0.88,r'Low $cAMP_{e}$ input', ha='center',va='center',
     transform = ax51.transAxes, color = 'k', fontsize= text_size)
ax51.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax51.axvspan(15, 30, alpha=0.2, color=simcolor)
ax51.set_xlim([0,30]); ax51.set_ylim([-0.4,1.6])
ax51.text(-0.06 , 1.2, 'F', ha='center',va='center',
         transform = ax51.transAxes, color = simcolor, fontsize=abcd_font_size)
ax51.set_title('IFFL',color = kaminocolor, size=title_font_size)

ax52= fig.add_subplot(grid[4, 2])
for count in range(5):
    ax52.plot(t_plot_Kamino,y_traces_norm[1,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax52.plot(t_plot_Kamino,y_traces_norm_mean[1,:], color= kaminocolor, linewidth=trace_width)
ax52.text(0.73,0.75,r'Intermediate $cAMP_{e}$'+'\n input', ha='center',va='center',
     transform = ax52.transAxes, color = 'k', fontsize= text_size)
ax52.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax52.axvspan(15, 30, alpha=0.2, color='g')
ax52.set_xlim([0,30]); ax52.set_ylim([-0.4,1.6])

ax53= fig.add_subplot(grid[5, 2])
for count in range(5):
    ax53.plot(t_plot_Kamino,y_traces_norm[2,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax53.plot(t_plot_Kamino,y_traces_norm_mean[2,:], color= kaminocolor, linewidth=trace_width)
ax53.text(0.73,0.88,r'High $cAMP_{e}$ input', ha='center',va='center',
     transform = ax53.transAxes, color = 'k', fontsize= text_size)
ax53.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
ax53.axvspan(15, 30, alpha=0.2, color=simcolor)
ax53.set_xlim([0,30]); ax53.set_ylim([-0.4,1.6]) 

fig.text(0.01, 0.36, r'$cAMP_{i}$', fontsize=label_font_size,va='center', rotation='vertical')
fig.text(0.5, 0, 'Time', fontsize=label_font_size, ha='center')
plt.subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.07, wspace=0.15, hspace=0.6)

# plt.subplots_adjust(left=0.1, right=0.96, top=0.93, bottom=0.04, wspace=0.4, hspace=0.72)
# grid.tight_layout(fig5,rect=[0, 0, 1, 1]) 
fig.savefig('fig7_pop_add_cAMP_tight_200530.png', bbox_inches='tight') 

#%% Fig 8 pop FR with SC noise
# # Experimental data
# PopRateExp = np.array([[0.129660088,0.171190476,0.152777778,0.171016069,0.132986111,0.149074074,0.144444444,0.177170868,0.162719633,0.178571429],
# [0.065434419,0.131937322,0.131658497,0.156811146,0.147619048,0.162799043,0.147368421,0.171820616,0.11875,0.1],
# [0.075416667,0.152597403,0.138703704,0.161497326,0.183333333,0.157142857,0.132373581,0.139815266,0.117989418,0.103125],
# [0.057236842,0.164138177,0.13287037,0.165909091,0.126666667,0.166827485,0.086222222,0.065789474,0.047584541,0.055555556],
# [0.019047619,0.103693161,0.122423896,0.100198413,0.106683375,0.091592443,0.034122807,0.042105263,0.031818182,0],
# [0.014285714,0.026428571,0.031578947,0.049938272,0.038095238,0.00625,0.010526316,0,0,0],
# [0,0.02375,0.022222222,0.050595238,0,0,0,0,0,0]])
# JExp = np.linspace(1,10, num=10) # np.array([0,1,2,4,6,8,10,15,16,20])
# RhoExp = np.linspace(1,7, num=7) #  np.array([0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125])
# np.savez('Gregor2010_pop_firing_rate.npz', PopRateExp = PopRateExp,
#          JExp = JExp, RhoExp = RhoExp)
  # load experimental data
npzfile = np.load('../exp_data/Gregor2010_pop_firing_rate.npz')
PopRateExp = npzfile['PopRateExp']
JExp = npzfile['JExp']; RhoExp = npzfile['RhoExp']

# Simulation outputs
# Goldbeter 1987 from parallel computing outputs
OUT_path = 'C:/Users/ellin/Documents/GitHub/dictymodels/model_outputs/'
Gold_OUT = np.load(OUT_path + 'pop_FR_Goldbeter_200311_hnorm_dt0.001_noise0ParamLen40.npz')
kc_arr_Gold =  Gold_OUT['kc_arr']
h_arr_Gold =  Gold_OUT['h_arr']
oneoverh_arr_Gold = 1/h_arr_Gold
pop_rate_Gold = Gold_OUT['pop_rate_Goldbeter']

Gold_OUT_noise10 = np.load(OUT_path + 'pop_FR_Goldbeter_200311_hnorm_dt0.001_noise10ParamLen40.npz')
kc_arr_Gold_noise10 =  Gold_OUT_noise10['kc_arr']
h_arr_Gold_noise10 =  Gold_OUT_noise10['h_arr']
oneoverh_arr_Gold_noise10 = 1/h_arr_Gold_noise10
pop_rate_Gold_noise10 = Gold_OUT_noise10['pop_rate_Goldbeter']
#  Maeda Loomis 2004 from parallel computing outputs
Maeda_OUT = np.load(OUT_path +'pop_FR_Maeda_200331_hnorm_dt0.0001_noise0ParamLen25p.npz')
rho_arr_Maeda =  Maeda_OUT['rho_arr']
gamma_arr_Maeda =  Maeda_OUT['gamma_arr']
pop_rate_Maeda = Maeda_OUT['pop_rate_Maeda']

Maeda_OUT_noise1 = np.load(OUT_path +'pop_FR_Maeda_200331_hnorm_dt0.0001_noise1ParamLen25p.npz')
rho_arr_Maeda_noise1 =  Maeda_OUT_noise1['rho_arr']
gamma_arr_Maeda_noise1 =  Maeda_OUT_noise1['gamma_arr']
pop_rate_Maeda_noise1 = Maeda_OUT_noise1['pop_rate_Maeda']

# Phase oscillator from parallel computing outputs
Gregor_OUT_noise = np.load(OUT_path +'pop_fire_rate_Gregor_OUT_191026_noise0.002.npz')
k_arr_Gregor_noise =  Gregor_OUT_noise['k_arr']
rho_arr_Gregor_noise =  Gregor_OUT_noise['rho_arr']
pop_rate_Gregor_noise = Gregor_OUT_noise['pop_rate_Gregor']

Gregor_OUT = np.load(OUT_path +'pop_fire_rate_Gregor_OUT_191027_noise0.npz')
k_arr_Gregor =  Gregor_OUT['k_arr']
rho_arr_Gregor =  Gregor_OUT['rho_arr']
pop_rate_Gregor = Gregor_OUT['pop_rate_Gregor']

# IPNFB
Sgro_no_noise_OUT = np.load(OUT_path +'pop_fire_rate_Sgro_OUT_200414_same_init_cond_ttot_25.0_dt0.005_sigma0_dir_cpl0.npz')
j_arr_Sgro_no_noise =  Sgro_no_noise_OUT['j_arr']
rho_arr_Sgro_no_noise =  Sgro_no_noise_OUT['rho_arr']
pop_rate_Sgro_no_noise =  Sgro_no_noise_OUT['pop_rate_Sgro']

Sgro_low_noise_OUT = np.load(OUT_path +'pop_fire_rate_Sgro_OUT_200414_same_init_cond_ttot_25.0_dt0.005_sigma0.1_dir_cpl0.npz')
j_arr_Sgro_low_noise =  Sgro_low_noise_OUT['j_arr']
rho_arr_Sgro_low_noise =  Sgro_low_noise_OUT['rho_arr']
pop_rate_Sgro_low_noise =  Sgro_low_noise_OUT['pop_rate_Sgro']

Sgro_regular_noise_OUT = np.load(OUT_path +'pop_fire_rate_Sgro_OUT_200413_same_init_cond_ttot_25.0_dt0.005_sigma0.15_dir_cpl0.npz')
j_arr_Sgro_regular_noise =  Sgro_regular_noise_OUT['j_arr']
rho_arr_Sgro_regular_noise =  Sgro_regular_noise_OUT['rho_arr']
pop_rate_Sgro_regular_noise =  Sgro_regular_noise_OUT['pop_rate_Sgro']

# IFFL from parallel computing outputs
Kamino_OUT = np.load(OUT_path +'pop_FR_Kamino_200401_dt0.001_noise_0_PrmLen_25PkFindThr0.03.npz')
rho_arr_Kamino =  Kamino_OUT['rho_arr']
gamma_arr_Kamino =  Kamino_OUT['gamma_arr']
pop_rate_Kamino = Kamino_OUT['pop_rate_Kamino']

Kamino_OUT_noise = np.load(OUT_path +'pop_FR_Kamino_200401_dt0.001_noise_0.01_PrmLen_25PkFindThr0.03.npz')
rho_arr_Kamino_noise =  Kamino_OUT_noise['rho_arr']
gamma_arr_Kamino_noise =  Kamino_OUT_noise['gamma_arr']
pop_rate_Kamino_noise = Kamino_OUT_noise['pop_rate_Kamino']

#%%
title_font_size = 18
label_font_size = 20
sublabel_font_size = 18
tick_font_size = 18
trace_width = 2
abcd_font_size = 28

fig3 = plt.figure(figsize=(cm2inch(48),cm2inch(24)))
grid = plt.GridSpec(3,5, wspace=0.35, hspace=0.6)

ax0= fig3.add_subplot(grid[0,0])
ax0.set_xticks([0,2,4,6,8]); 
ax0.set_xticklabels([1,3,5,7,9],fontsize=tick_font_size)
# ax0.set_yticks([0,1,2,3,4,5,6,7]); 
# ax0.set_yticklabels(['1/2','1/4','1/8','1/16','1/32','1/64','1/128'],fontsize=tick_font_size-3)
ax0.set_yticks([0,2,4,6]); 
ax0.set_yticklabels(['1/2','1/8','1/32','1/128'],fontsize=tick_font_size)

ax0.set_title('Experiment', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax0.set_xlabel('Flow Rate (mL/min)', size=sublabel_font_size)
ax0.set_ylabel('Cell Density(mML)', size=sublabel_font_size)
heatmap = ax0.imshow(PopRateExp, cmap='jet') # cmap='jet'
x=[3.5,4.5,5.5,7.5,9.5]
[ax0.axvline(_x, color='white',linewidth=trace_width) for _x in x]
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax0,ticks=[0,0.05,0.1,0.15]);
cbar.ax.tick_params(labelsize = tick_font_size) 
cbar.set_label( 'cAMP pulses/min',size=tick_font_size)
#ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.text(-0.6 , 1.6, 'A',
         horizontalalignment='center',verticalalignment='center',
         transform = ax0.transAxes, color = 'k', fontsize=abcd_font_size)

ax1= fig3.add_subplot(grid[1,0],xticks=[0,25,50,75,100])
heatmap = ax1.pcolor(kc_arr_Gold_noise10, oneoverh_arr_Gold_noise10, pop_rate_Gold_noise10.transpose(), cmap='jet') # cmap='jet'
# ax1.set_xscale('log');
ax1.set_yscale('log')
heatmap.set_clim(0,1.5)
cbar=fig3.colorbar(heatmap, ax=ax1, ticks=[0,0.5,1,1.5])
cbar.ax.tick_params(labelsize = tick_font_size) 
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Receptor Desensitization\nwith noise', color = goldcolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax1.text(-0.6 , 1.2, 'B',
         horizontalalignment='center',verticalalignment='center',
         transform = ax1.transAxes, color = 'k', fontsize=abcd_font_size)

ax1lower= fig3.add_subplot(grid[2,0],xticks=[0,25,50,75,100])
heatmap = ax1lower.pcolor(kc_arr_Gold, oneoverh_arr_Gold, pop_rate_Gold.transpose(), cmap='jet') # cmap='jet'
# ax1.set_xscale('log');
ax1lower.set_yscale('log')
heatmap.set_clim(0,1.5)
cbar=fig3.colorbar(heatmap, ax=ax1lower, ticks=[0,0.5,1,1.5])
cbar.ax.tick_params(labelsize = tick_font_size) 
ax1lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1lower.set_title('Receptor Desensitization\nw/o noise', color = goldcolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax2= fig3.add_subplot(grid[1,1],xticks=[0,25,50,75,100])
heatmap = ax2.pcolor(gamma_arr_Maeda_noise1, rho_arr_Maeda_noise1, pop_rate_Maeda_noise1.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.55)
ax2.set_yscale('log')
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('CDINFB\nwith noise', color = maedacolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax2.text(-0.3 , 1.2, 'C',
         horizontalalignment='center',verticalalignment='center',
         transform = ax2.transAxes, color = 'k', fontsize=abcd_font_size)

ax2lower= fig3.add_subplot(grid[2,1],xticks=[0,25,50,75,100])
heatmap = ax2lower.pcolor(gamma_arr_Maeda, rho_arr_Maeda, pop_rate_Maeda.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.55)
ax2lower.set_yscale('log')
cbar=fig3.colorbar(heatmap, ax=ax2lower);cbar.ax.tick_params(labelsize = tick_font_size) 
ax2lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2lower.set_title('CDINFB\nw/o noise', color= maedacolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax3= fig3.add_subplot(grid[1,2], xticks=[0,25,50,75,100])
heatmap = ax3.pcolor(k_arr_Gregor_noise, rho_arr_Gregor_noise, pop_rate_Gregor_noise.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,1.2)
ax3.set_yscale('log')
cbar=fig3.colorbar(heatmap, ax=ax3);cbar.ax.tick_params(labelsize = tick_font_size) 
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Phase Oscillator\nwith noise',color= gregorcolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax3.text(-0.3 , 1.2, 'D',
         horizontalalignment='center',verticalalignment='center',
         transform = ax3.transAxes, color = 'k', fontsize=abcd_font_size)

ax3lower = fig3.add_subplot(grid[2,2], xticks=[0,25,50,75,100])
heatmap = ax3lower.pcolor(k_arr_Gregor, rho_arr_Gregor, pop_rate_Gregor.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,1.2)
ax3lower.set_yscale('log')
cbar=fig3.colorbar(heatmap, ax=ax3lower);cbar.ax.tick_params(labelsize = tick_font_size) 
ax3lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3lower.set_title('Phase Oscillator\nw/o noise',color= gregorcolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# Sgro regular noise (sig = 0.15)
ax4= fig3.add_subplot(grid[1,3],xticks=[0,0.5,1])
heatmap = ax4.pcolor(j_arr_Sgro_regular_noise, rho_arr_Sgro_regular_noise, 
                     pop_rate_Sgro_regular_noise.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.55)
ax4.set_yscale('log'); # ax4.set_ylim([10**(-5),10**(-3)]); 
cbar=fig3.colorbar(heatmap, ax=ax4);cbar.ax.tick_params(labelsize = tick_font_size) 
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('IPNFB\nwith noise', color= sgrocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax4.text(-0.3 , 1.2, 'E',
         horizontalalignment='center',verticalalignment='center',
         transform = ax4.transAxes, color = 'k', fontsize=abcd_font_size)
# ax4.set_xticklabels([0,0.25,0.5,0.75,1], rotation=45,fontsize=tick_font_size)

# Sgro w/o noise  (sig = 0.0)
ax4lower= fig3.add_subplot(grid[2,3],xticks=[0,0.5,1])
heatmap = ax4lower.pcolor(j_arr_Sgro_no_noise, rho_arr_Sgro_no_noise, 
                          pop_rate_Sgro_no_noise.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.55)
ax4lower.set_yscale('log'); # ax4lower.set_ylim([10**(-5),10**(-3)]); 
cbar=fig3.colorbar(heatmap, ax=ax4lower);cbar.ax.tick_params(labelsize = tick_font_size) 
ax4lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4lower.set_title('IPNFB\nw/o noise', color= sgrocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# # Sgro low noise  (sig = 0.1)
# ax4lower= fig3.add_subplot(grid[2,3],xticks=[0,0.25,0.5,0.75,1])
# heatmap = ax4lower.pcolor(j_arr_Sgro_low_noise, rho_arr_Sgro_low_noise, 
#                           pop_rate_Sgro_low_noise.transpose(), cmap='jet') # cmap='jet'
# heatmap.set_clim(0,0.6)
# ax4lower.set_yscale('log'); ax4lower.set_ylim([10**(-5),10**(-3)]); 
# cbar=fig3.colorbar(heatmap, ax=ax4lower);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax4lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
# ax4lower.set_title('IPNFB\nlow noise', color=mycolors[5],fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})


ax5= fig3.add_subplot(grid[1,4], xticks=[0,25,50,75,100])
heatmap = ax5.pcolor(gamma_arr_Kamino_noise, rho_arr_Kamino_noise, 
                     pop_rate_Kamino_noise.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.65)
ax5.set_yscale('log')
cbar=fig3.colorbar(heatmap, ax=ax5);cbar.ax.tick_params(labelsize = tick_font_size) 
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5.set_title('IFFL\nwith noise', color= kaminocolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax5.text(-0.3 , 1.2, 'F',
         horizontalalignment='center',verticalalignment='center',
         transform = ax5.transAxes, color = 'k', fontsize=abcd_font_size)

# insets
#ax7= fig3.add_axes([0.77,0.13,0.08,0.16])
ax5in= fig3.add_axes([0.835,0.41,0.042,0.082])
heatmap = ax5in.pcolor(gamma_arr_Kamino_noise, rho_arr_Kamino_noise,
                       pop_rate_Kamino_noise.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.65)
ax5in.set_yscale('log');ax5in.set_xscale('log');
ax5in.set_xticks([]) ; ax5in.set_yticks([]) 
ax5in.spines['bottom'].set_color('white');ax5in.spines['top'].set_color('white')
ax5in.spines['left'].set_color('white');ax5in.spines['right'].set_color('white')

ax5lower= fig3.add_subplot(grid[2,4], xticks=[0,25,50,75,100])
heatmap = ax5lower.pcolor(gamma_arr_Kamino, rho_arr_Kamino, 
                          pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.65)
ax5lower.set_yscale('log')
cbar=fig3.colorbar(heatmap, ax=ax5lower);cbar.ax.tick_params(labelsize = tick_font_size) 
ax5lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5lower.set_title('IFFL\nw/o noise', color = kaminocolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax5lowerin= fig3.add_axes([0.835,0.115,0.042,0.082])
heatmap = ax5lowerin.pcolor(gamma_arr_Kamino, rho_arr_Kamino,
                       pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.65)
ax5lowerin.set_yscale('log'); ax5lowerin.set_xscale('log');
ax5lowerin.set_xticks([]) ; ax5lowerin.set_yticks([]) 
ax5lowerin.spines['bottom'].set_color('white');ax5lowerin.spines['top'].set_color('white')
ax5lowerin.spines['left'].set_color('white');ax5lowerin.spines['right'].set_color('white')


fig3.text(0.51, 0.045, 'Dilution Rate, A.U.',fontsize=label_font_size, ha='center')
fig3.text(0.08, 0.35, 'Population Density, A.U.',fontsize=label_font_size, va='center', rotation='vertical')

# plt.subplots_adjust(left=0.1, right=0.96, top=0.93, bottom=0.04, wspace=0.4, hspace=0.72)
# grid.tight_layout(fig5,rect=[0, 0, 1, 1]) 
#fig3.savefig('fig8_pop_firing_rate_tight_200530.png', bbox_inches='tight') 

#%% fig S9-12 population firing rate
# import population firing rate heatmap outputs
# Goldbeter 1987s
OUT_path = 'C:/Users/ellin/Documents/GitHub/dictymodels/model_outputs/'
Gold_OUT = np.load(OUT_path + 'pop_FR_Goldbeter_200311_hnorm_dt0.001_noise0ParamLen40.npz')
Gold_OUT_noise10 = np.load(OUT_path + 'pop_FR_Goldbeter_200311_hnorm_dt0.001_noise10ParamLen40.npz')
kc_arr_Gold_noise10 =  Gold_OUT_noise10['kc_arr']
h_arr_Gold_noise10 =  Gold_OUT_noise10['h_arr']
oneoverh_arr_Gold_noise10 = 1/h_arr_Gold_noise10
pop_rate_Gold_noise10 = Gold_OUT_noise10['pop_rate_Goldbeter']
#  Maeda Loomis 2004
Maeda_OUT_noise1 = np.load(OUT_path +'pop_FR_Maeda_200331_hnorm_dt0.0001_noise1ParamLen25p.npz')
rho_arr_Maeda_noise1 =  Maeda_OUT_noise1['rho_arr']
gamma_arr_Maeda_noise1 =  Maeda_OUT_noise1['gamma_arr']
pop_rate_Maeda_noise1 = Maeda_OUT_noise1['pop_rate_Maeda']
# Phase oscillator 
Gregor_OUT_noise = np.load(OUT_path +'pop_fire_rate_Gregor_OUT_191026_noise0.002.npz')
k_arr_Gregor_noise =  Gregor_OUT_noise['k_arr']
rho_arr_Gregor_noise =  Gregor_OUT_noise['rho_arr']
pop_rate_Gregor_noise = Gregor_OUT_noise['pop_rate_Gregor']
# IPNFB
Sgro_regular_noise_OUT = np.load(OUT_path +'pop_fire_rate_Sgro_OUT_200413_same_init_cond_ttot_25.0_dt0.005_sigma0.15_dir_cpl0.npz')
j_arr_Sgro_regular_noise =  Sgro_regular_noise_OUT['j_arr']
rho_arr_Sgro_regular_noise =  Sgro_regular_noise_OUT['rho_arr']
pop_rate_Sgro_regular_noise =  Sgro_regular_noise_OUT['pop_rate_Sgro']
# IFFL 
Kamino_OUT_noise = np.load(OUT_path +'pop_FR_Kamino_200401_dt0.001_noise_0.01_PrmLen_25PkFindThr0.03.npz')
rho_arr_Kamino_noise =  Kamino_OUT_noise['rho_arr']
gamma_arr_Kamino_noise =  Kamino_OUT_noise['gamma_arr']
pop_rate_Kamino_noise = Kamino_OUT_noise['pop_rate_Kamino']

# load parameter pair  npz output file 
import_npz(OUT_path+'pop_FR_supp_200426.npz',globals())
import_npz(OUT_path+'pop_FR_supp_200430_cont.npz',globals())
#%% Specify plottig sizes
title_font_size = 20
label_font_size = 20
sublabel_font_size = 18
tick_font_size = 16
text_size = 14
trace_width = 2
abcd_font_size = 28
m_size = 8; m_color = 'darkgrey'

SC_traces_idx = [0,2,4,6,8,10]

#%% Fig S7 Goldbeter
fig, ax = plt.subplots(nrows=2, ncols=4, figsize=(cm2inch(43),cm2inch(20)))
plt.subplots_adjust(left=0.1, right=0.97, top=0.9, bottom=0.1, 
                    wspace=0.24, hspace=0.6)
ax[0,0].axis('off'); ax[1,0].axis('off')

ax0 = fig.add_axes([0.06, 0.55, 0.17,0.33],xticks=[25,50,75,100])
heatmap = ax0.pcolor(kc_arr_Gold_noise10, oneoverh_arr_Gold_noise10, pop_rate_Gold_noise10.transpose(), cmap='jet') # cmap='jet'
# ax1.set_xscale('log');
ax0.set_yscale('log')
heatmap.set_clim(0,1.5)
ax0.set_xlabel('Dilution Rate, A.U.', size=sublabel_font_size-2)
ax0.set_ylabel('Cell Density, A.U.', size=sublabel_font_size-2)
cbar=fig.colorbar(heatmap, ax=ax0, ticks=[0,0.5,1,1.5])
cbar.ax.tick_params(labelsize = tick_font_size) 
cbar.set_label( 'population firing rate',size=tick_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_title('Receptor desensitization\nwith noise', color = goldcolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax0.plot(kc_arr[0:3],1/h_arr[0:3],'X',markersize =  m_size,color = m_color)
ax0.plot(kc_arr[3:],1/h_arr[3:],'v',markersize =  m_size,color = m_color)

for temp in range(6):
    row = int(temp/3); col = temp%3 + 1
    title = 'Density='+'{:#.3n}'.format(np.float64(1/h_arr[temp]))+\
            ', dilution rate='+ '{:#.3n}'.format(np.float64(kc_arr[temp]))+\
                '\nFR= '+'{:#.3n}'.format(np.float64(FR_gold[temp])) 
    for this_idx in SC_traces_idx:
        this_trace = b_traces_later[temp,this_idx,:]# /np.amax(SC_traces[this_idx,:])
        ax[row,col].plot(t_plot_Goldbeter_short,this_trace, 
                         color='darkgrey', alpha=0.6, linewidth=3)     
    # Plot population mean
    ax[row,col].plot(t_plot_Goldbeter_short, b_traces_mean_later[temp,:], 
                     color=goldcolor ,linewidth=2, label = 'cAMPi')
    ax[row,col].set_title(title, fontdict={'fontsize':sublabel_font_size})
    # ax1[row,col].set_xlabel('Time, A.U.', fontsize=sublabel_font_size); 
    # ax1.set_ylabel(r'cAMPi', color = 'k', fontsize=20 )
    ax[row,col].tick_params(axis = 'y', labelcolor = 'k')
    ax[row,col].tick_params( grid_linewidth = tick_font_size, labelsize = tick_font_size)
    ax[row,col].set_xlim([15,25])
ax0.text(-0.3, 1.21, 'A', ha='center',va='center',
         transform = ax0.transAxes, color = 'k', fontsize=abcd_font_size)
ax[0,1].text(-0.25, 1.15, 'B', ha='center',va='center',
         transform = ax[0,1].transAxes, color = 'k', fontsize=abcd_font_size)
fig.text(0.62, 0.02, ' Time, A.U.', ha='center', va='center',fontsize=label_font_size)
fig.text(0.29, 0.55, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)      
 
fig.savefig('figS7_popFR_Goldbeter_tight_200530.png', bbox_inches='tight') 

#%% Fig S8 Maeda
fig, ax = plt.subplots(nrows=2, ncols=4, figsize=(cm2inch(43),cm2inch(20)))
plt.subplots_adjust(left=0.1, right=0.97, top=0.9, bottom=0.1, 
                    wspace=0.24, hspace=0.6)
ax[0,0].axis('off'); ax[1,0].axis('off')

ax0 = fig.add_axes([0.06, 0.55, 0.17,0.33],xticks=[25,50,75,100])
heatmap = ax0.pcolor(gamma_arr_Maeda_noise1, rho_arr_Maeda_noise1, 
                     pop_rate_Maeda_noise1.transpose(), cmap='jet') # cmap='jet'
# ax1.set_xscale('log');
ax0.set_yscale('log')
heatmap.set_clim(0,0.55)
ax0.set_xlabel('Dilution Rate, A.U.', size=sublabel_font_size-2)
ax0.set_ylabel('Cell Density, A.U.', size=sublabel_font_size-2)
cbar=fig.colorbar(heatmap, ax=ax0, ticks=[0,0.5,1,1.5])
cbar.ax.tick_params(labelsize = tick_font_size) 
cbar.set_label( 'population firing rate',size=tick_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_title('CDINFB\nwith noise', color = maedacolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax0.plot(gamma_arr_Maeda[0:3],rho_arr_Maeda[0:3],'X',markersize =  m_size,color = m_color)
ax0.plot(gamma_arr_Maeda[3:],rho_arr_Maeda[3:],'v',markersize =  m_size,color = m_color)

for temp in range(6):
    row = int(temp/3); col = temp%3 + 1
    title = 'Density='+'{:#.3n}'.format(np.float64(rho_arr_Maeda[temp]))+\
            ', dilution rate='+ '{:#.3n}'.format(np.float64(gamma_arr_Maeda[temp]))+\
                '\nFR= '+'{:#.3n}'.format(np.float64(FR_Maeda[temp])) 
    for this_idx in SC_traces_idx:
        this_trace = cAMPi_traces_later [temp,this_idx,:]# /np.amax(SC_traces[this_idx,:])
        ax[row,col].plot(t_plot_Maeda_short,this_trace, 
                         color='darkgrey', alpha=0.6, linewidth=3)     
    # Plot population mean
    ax[row,col].plot(t_plot_Maeda_short, cAMPi_traces_mean_later [temp,:], 
                     color=maedacolor ,linewidth=2, label = 'cAMPi')
    ax[row,col].set_title(title, fontdict={'fontsize':sublabel_font_size})
    # ax1[row,col].set_xlabel('Time, A.U.', fontsize=sublabel_font_size); 
    # ax1.set_ylabel(r'cAMPi', color = 'k', fontsize=20 )
    ax[row,col].tick_params(axis = 'y', labelcolor = 'k')
    ax[row,col].tick_params( grid_linewidth = tick_font_size, labelsize = tick_font_size)
    ax[row,col].set_xlim([30,50])
ax0.text(-0.3, 1.21, 'A', ha='center',va='center',
         transform = ax0.transAxes, color = 'k', fontsize=abcd_font_size)
ax[0,1].text(-0.25, 1.15, 'B', ha='center',va='center',
         transform = ax[0,1].transAxes, color = 'k', fontsize=abcd_font_size)
fig.text(0.62, 0.02, ' Time, A.U.', ha='center', va='center',fontsize=label_font_size)
fig.text(0.29, 0.55, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)       

fig.savefig('figS8_popFR_Maeda_tight_200530.png', bbox_inches='tight') 

#%% Fig S9 Phase oscillator
fig, ax = plt.subplots(nrows=2, ncols=4, figsize=(cm2inch(43),cm2inch(20)))
plt.subplots_adjust(left=0.1, right=0.97, top=0.9, bottom=0.1, 
                    wspace=0.24, hspace=0.6)
ax[0,0].axis('off'); ax[1,0].axis('off')

ax0 = fig.add_axes([0.06, 0.55, 0.17,0.33],xticks=[25,50])
heatmap = ax0.pcolor(k_arr_Gregor_noise, rho_arr_Gregor_noise, 
                     pop_rate_Gregor_noise.transpose(), cmap='jet') # cmap='jet'
# ax1.set_xscale('log');
ax0.set_xlim([0,50])
ax0.set_yscale('log')
heatmap.set_clim(0,1.2)
ax0.set_xlabel('Dilution Rate, A.U.', size=sublabel_font_size-2)
ax0.set_ylabel('Cell Density, A.U.', size=sublabel_font_size-2)
cbar.ax.tick_params(labelsize = tick_font_size) 
cbar.set_label( 'population firing rate',size=tick_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_title('Greor 2010\nwith noise', color = gregorcolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax0.plot(k_arr_Gregor[0:3], rho_arr_Gregor[0:3],'X',markersize =  m_size,color = m_color)
ax0.plot(k_arr_Gregor[3:], rho_arr_Gregor[3:],'v',markersize =  m_size,color = m_color)

for temp in range(6):
    row = int(temp/3); col = temp%3 + 1
    title = 'Density='+'{:#.3n}'.format(np.float64(rho_arr_Gregor[temp]))+\
            ', dilution rate='+ '{:#.3n}'.format(np.float64(k_arr_Gregor[temp]))+\
                '\nFR= '+'{:#.3n}'.format(np.float64(FR_Gregor[temp])) 
    for this_idx in SC_traces_idx:
        this_trace = cAMPi_Gregor_later [temp,this_idx,:]# /np.amax(SC_traces[this_idx,:])
        ax[row,col].plot(t_plot_Gregor_later,this_trace, 
                         color='darkgrey', alpha=0.6, linewidth=3)     
    # Plot population mean
    ax[row,col].plot(t_plot_Gregor_later, cAMPi_Gregor_mean_later [temp,:], 
                     color=gregorcolor ,linewidth=2, label = 'cAMPi')
    ax[row,col].set_title(title, fontdict={'fontsize':sublabel_font_size})
    # ax1[row,col].set_xlabel('Time, A.U.', fontsize=sublabel_font_size); 
    # ax1.set_ylabel(r'cAMPi', color = 'k', fontsize=20 )
    ax[row,col].tick_params(axis = 'y', labelcolor = 'k')
    ax[row,col].tick_params( grid_linewidth = tick_font_size, labelsize = tick_font_size)
    ax[row,col].set_xlim([15,25])
ax0.text(-0.3, 1.21, 'A', ha='center',va='center',
         transform = ax0.transAxes, color = 'k', fontsize=abcd_font_size)
ax[0,1].text(-0.25, 1.15, 'B', ha='center',va='center',
         transform = ax[0,1].transAxes, color = 'k', fontsize=abcd_font_size)
fig.text(0.62, 0.02, ' Time, A.U.', ha='center', va='center',fontsize=label_font_size)
fig.text(0.29, 0.55, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)       
fig.savefig('figS9_popFR_Gregor_tight_200530.png', bbox_inches='tight') 

#%% Fig S10 IPNFB
m_size = 8; m_color = 'w'

fig, ax = plt.subplots(nrows=2, ncols=4, figsize=(cm2inch(43),cm2inch(20)))
plt.subplots_adjust(left=0.1, right=0.97, top=0.9, bottom=0.1, 
                    wspace=0.24, hspace=0.6)
ax[0,0].axis('off'); ax[1,0].axis('off')

ax0 = fig.add_axes([0.06, 0.55, 0.17,0.33],xticks=[0,0.5,1])
heatmap = ax0.pcolor(j_arr_Sgro_regular_noise, rho_arr_Sgro_regular_noise, 
                     pop_rate_Sgro_regular_noise.transpose(), cmap='jet') # cmap='jet'
# ax1.set_xscale('log');
ax0.set_yscale('log')
heatmap.set_clim(0,0.55)
ax0.set_xlabel('Dilution Rate, A.U.', size=sublabel_font_size-2)
ax0.set_ylabel('Cell Density, A.U.', size=sublabel_font_size-2)
cbar.ax.tick_params(labelsize = tick_font_size) 
cbar.set_label( 'population firing rate',size=tick_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_title('IPNFB\nwith noise', color = sgrocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax0.plot(j_arr_Sgro[0:3], rho_arr_Sgro[0:3],'X',markersize =  m_size,color = m_color)
ax0.plot(j_arr_Sgro[3:], rho_arr_Sgro[3:],'v',markersize =  m_size,color = m_color)

for temp in range(6):
    row = int(temp/3); col = temp%3 + 1
    title = 'Density='+'{:#.3n}'.format(np.float64(rho_arr_Sgro[temp]))+\
            ',\ndilution rate='+ '{:#.3n}'.format(np.float64(j_arr_Sgro[temp]))+\
                ', FR= '+'{:#.3n}'.format(np.float64(FR_Sgro[temp])) 
    for this_idx in SC_traces_idx:
        this_trace = A_traces_later[temp,this_idx,:]# /np.amax(SC_traces[this_idx,:])
        ax[row,col].plot(t_plot_Sgro_later,this_trace, 
                         color='darkgrey', alpha=0.6, linewidth=3)     
    # Plot population mean
    ax[row,col].plot(t_plot_Sgro_later, A_traces_mean_later[temp,:], 
                     color=sgrocolor ,linewidth=2, label = 'cAMPi')
    ax[row,col].set_title(title, fontdict={'fontsize':sublabel_font_size})
    ax[row,col].tick_params(axis = 'y', labelcolor = 'k')
    ax[row,col].tick_params( grid_linewidth = tick_font_size, labelsize = tick_font_size)
    ax[row,col].set_xlim([t_plot_Sgro_later[0],t_plot_Sgro_later[-1]])
ax0.text(-0.3, 1.21, 'A', ha='center',va='center',
         transform = ax0.transAxes, color = 'k', fontsize=abcd_font_size)
ax[0,1].text(-0.25, 1.15, 'B', ha='center',va='center',
         transform = ax[0,1].transAxes, color = 'k', fontsize=abcd_font_size)
fig.text(0.62, 0.02, ' Time, A.U.', ha='center', va='center',fontsize=label_font_size)
fig.text(0.29, 0.55, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)       

fig.savefig('figS10_popFR_Sgro_tight_200530.png', bbox_inches='tight') 

#%% Fig S11 IFFL
m_size = 8; m_color = 'lightgrey'

fig, ax = plt.subplots(nrows=2, ncols=4, figsize=(cm2inch(45),cm2inch(20)))
plt.subplots_adjust(left=0.1, right=0.97, top=0.9, bottom=0.1, 
                    wspace=0.26, hspace=0.6)
ax[0,0].axis('off'); ax[1,0].axis('off')

ax0 = fig.add_axes([0.06, 0.55, 0.17,0.33],xticks=[0,25,50,75,100])
heatmap = ax0.pcolor(gamma_arr_Kamino_noise, rho_arr_Kamino_noise, 
                     pop_rate_Kamino_noise.transpose(), cmap='jet') # cmap='jet'
ax0.set_xscale('log');
ax0.set_yscale('log')
heatmap.set_clim(0,0.65)
ax0.set_xlabel('Dilution Rate, A.U.', size=sublabel_font_size-2)
ax0.set_ylabel('Cell Density, A.U.', size=sublabel_font_size-2)
cbar.ax.tick_params(labelsize = tick_font_size) 
cbar.set_label( 'population firing rate',size=tick_font_size)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_title('IFFL\nwith noise', color = kaminocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax0.plot(gamma_arr_Kamino[0:3], rho_arr_Kamino[0:3],'X',markersize =  m_size,color = m_color)
ax0.plot(gamma_arr_Kamino[3:], rho_arr_Kamino[3:],'v',markersize =  m_size,color = m_color)

for temp in range(6):
    row = int(temp/3); col = temp%3 + 1
    title = 'Density='+'{:#.3n}'.format(np.float64(rho_arr_Kamino[temp]))+\
            ',\ndilution rate='+ '{:#.3n}'.format(np.float64(gamma_arr_Kamino[temp]))+\
                ', FR= '+'{:#.3n}'.format(np.float64(FR_Kamino[temp])) 
    for this_idx in SC_traces_idx:
        this_trace = y_traces_later[temp,this_idx,:]# /np.amax(SC_traces[this_idx,:])
        ax[row,col].plot(t_plot_Kamino_short,this_trace, 
                         color='darkgrey', alpha=0.6, linewidth=3)     
    # Plot population mean
    ax[row,col].plot(t_plot_Kamino_short, y_traces_mean_later[temp,:], 
                     color=kaminocolor ,linewidth=2, label = 'cAMPi')
    ax[row,col].set_title(title, fontdict={'fontsize':sublabel_font_size})
    ax[row,col].tick_params(axis = 'y', labelcolor = 'k')
    ax[row,col].tick_params( grid_linewidth = tick_font_size, labelsize = tick_font_size)
    ax[row,col].set_xlim([20,40])
ax0.text(-0.3, 1.21, 'A', ha='center',va='center',
         transform = ax0.transAxes, color = 'k', fontsize=abcd_font_size)
ax[0,1].text(-0.25, 1.15, 'B', ha='center',va='center',
         transform = ax[0,1].transAxes, color = 'k', fontsize=abcd_font_size)
fig.text(0.62, 0.02, ' Time, A.U.', ha='center', va='center',fontsize=label_font_size)
fig.text(0.29, 0.55, r'$cAMP_{i}$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)       

fig.savefig('figS11_popFR_Kamino_tight_200530.png', bbox_inches='tight') 

#%% Fig S13 Phase oscillator pop oscillations with and w/o noise
import_npz('../model_outputs/pop_osc_Gregor_noise0_200604.npz',globals())
import_npz('../model_outputs/pop_osc_Gregor_200604.npz',globals())

SC_traces_idx = [1,5,10,15,20,25,30]
label_font_size=20; trace_width=3; tick_font_size=18
title_font_size = 22

fig3 = plt.figure(figsize=(8,6))
grid = plt.GridSpec(2,1, wspace=0.15, hspace=0.4)
ax1=fig3.add_subplot(grid[0, 0])
ax1.plot(t_plot_Gregor,gregor_campCyto_trace_mean_norm,\
         linewidth=trace_width, color = 'g')
for this_idx in SC_traces_idx:
    this_trace = gregor_campCyto_trace_norm[this_idx,:]
    ax1.plot(t_plot_Gregor,this_trace, color='grey',alpha=0.6, linewidth=2)  
# ax1.set_xlabel('Time')
# ax1.set_ylabel(r'$cAMP_{i}$')
ax1.set_xlim([0.5*np.amax(t_plot_Gregor), np.amax(t_plot_Gregor)])
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title(r'With noise',fontdict={'fontsize': label_font_size})

ax2=fig3.add_subplot(grid[1, 0])
ax2.plot(t_plot_Gregor,gregor_campCyto_trace_mean_norm_noise0,\
         linewidth=trace_width, color = 'g')
for this_idx in SC_traces_idx:
    this_trace = gregor_campCyto_trace_norm_noise0[this_idx,:]
    ax2.plot(t_plot_Gregor,this_trace, color='grey',alpha=0.6, linewidth=2)  
ax2.set_xlim([0.5*np.amax(t_plot_Gregor), np.amax(t_plot_Gregor)])
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title(r'Without noise',fontdict={'fontsize': label_font_size})

fig3.text(0.05, 0.5,r'$cAMP_{i}$, A.U.', ha='center', va='center',rotation='vertical',
          fontsize=title_font_size)
fig3.text(0.5, 0.01, 'Time, A.U.', ha='center', 
          va='center',fontsize= title_font_size)
fig3.text(0.5, 1,'Phase oscillator population oscillations',ha='center', va='center',
          color= 'k',fontsize=title_font_size)
plt.show()

fig3.savefig('Fig13_200604.png', bbox_inches='tight')#
