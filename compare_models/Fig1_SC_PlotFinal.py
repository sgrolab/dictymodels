# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 11:39:13 2020

@author: cqhuyan
plotting script for figure 2- single cell behaviors

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
# original version colors
# colors = np.array([[69,183,204],
#                   [255,162,111],
#                   [10,68,12],
#                   [110,60,178],
#                   [222,76,105],
#                   [59,34,255], # color for experiments
#                   [35,145,40]])# color for models
# goldcolor = (colors[0,:]+(255-colors[0,:])*0.1)/255 # tint factor, larger means more "whiteness" to the original color
# maedacolor = (colors[1,:]+(255-colors[1,:])*0)/255
# gregorcolor = (colors[2,:]+(255-colors[2,:])*0.4)/255
# sgrocolor = (colors[3,:]+(255-colors[3,:])*0.3)/255
# kaminocolor = (colors[4,:]+(255-colors[4,:])*0.2)/255
# colorcombo1
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

#%% single cell and oscillations

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
ABCD_font_size = 24
abcd_font_size = 20
label_font_size=16
title_font_size = 16
sublabel_font_size = 14
trace_width=2
tick_font_size=14

fig = plt.figure(figsize=(18,24))
grid = plt.GridSpec(11, 11, wspace=0.8, hspace=0.3)
# grid = plt.GridSpec(9, 2, wspace=0.3, hspace=1.8)

axA01= fig.add_subplot(grid[0, 0:3])
axA01.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 1 FRET Trace (1nM)"],
                              linewidth=trace_width,color='k')
axA01.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 2 FRET Trace (1nM)"],
                               linewidth=trace_width,color='dimgrey')
axA01.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 3 FRET Trace (1nM)"],
                               linewidth=trace_width,color='darkgrey')
# axA01.set_ylabel(r'FRET Signal, A.U.',fontsize=label_font_size)
axA01.set_xlabel('Time (min)',fontsize=label_font_size)
axA01.axvline(x=5, ls='--', linewidth=trace_width, color = expcolor) #dashed line at 5 (cAMP onset)
axA01.set_ylim([-0.1,0.7]); axA01.set_xlim([0, 30])
axA01.tick_params(axis='both', which='major', labelsize=tick_font_size)
axA01.set_title('Experiment',color = 'k', fontsize = title_font_size)
axA01.text(0.7,0.8,' 1nM cAMP', ha='center',va='center',
         transform = axA01.transAxes, color = 'k', fontsize=label_font_size)
axA01.text(-0.06, 1.4, '(a)',
         ha='center',va='center',
         transform = axA01.transAxes, color = expcolor, fontsize=abcd_font_size)
axA01.text(-0.2, 2, 'A',
         ha='center',va='center',
         transform = axA01.transAxes, color = 'k', fontsize=ABCD_font_size)

axA02= fig.add_subplot(grid[2:4, 0:3])
#axA2.axvline(x = 5 , ymin=-0.2, ymax = 1.2, ls = '--', 
#            linewidth=trace_width, color=mycolors[6])
axA02.plot(t_plot_Goldbeter, b_trace,linewidth=trace_width, color = goldcolor,
         label='Receptor desensitization')
axA02.plot(t_plot_Maeda, cAMPi_trace,linewidth=trace_width,color = maedacolor,
         label='Coupled direct and indirect feedback')
axA02.tick_params(axis='both', which='major', labelsize=tick_font_size)
axA02.plot(t_plot_Gregor,gregor_campCyto_trace,linewidth=trace_width, 
         color= gregorcolor,label='Phase oscillator')
axA02.plot(t_plot_Sgro, A_trace_plot,linewidth=trace_width,
         color=sgrocolor, label='Interlocking positive-negative feedback')
axA02.plot(t_plot_Kamino, y_trace,linewidth=trace_width, 
         ls = '--',color=kaminocolor, label='IFFL')
axA02.set_title('Simulation',color = 'k', fontsize = title_font_size)
axA02.text(-0.06, 1.15, '(c)',
         ha='center',va='center',
         transform = axA02.transAxes, color = simcolor, fontsize=abcd_font_size)
axA02.set_ylim([-0.4,1.2])#([-0.2,1.2]) # adaptive spikes
axA02.set_xlim([0,6])
# axA02.set_xlabel('Time, A.U.',fontsize=label_font_size)
axA02.tick_params(axis='both', which='major', labelsize=tick_font_size)

# leg = axA2.legend()
axA02.legend(bbox_to_anchor=(2.15,-0.3),ncol=3,prop={'size': 11})


axA1= fig.add_subplot(grid[0, 3:6])
axA1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 1 FRET Trace (10uM)"],
                              linewidth=trace_width,color='k')
axA1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 2 FRET Trace (10uM)"],
                               linewidth=trace_width,color='dimgrey')
axA1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 3 FRET Trace (10uM)"],
                               linewidth=trace_width,color='darkgrey')

axA1.set_xlabel('Time (min)',fontsize=label_font_size)
axA1.axvline(x=5, ls='--', linewidth=trace_width, color = expcolor) #dashed line at 5 (cAMP onset)
axA1.set_ylim([-0.1,0.7]); axA1.set_xlim([0, 30])
axA1.tick_params(axis='both', which='major', labelsize=tick_font_size)
axA1.set_title('Experiment',color = 'k', fontsize = title_font_size)
axA1.text(0.7,0.8,' 10uM cAMP', ha='center',va='center',
         transform = axA1.transAxes, color = 'k', fontsize=label_font_size)
axA1.text(-0.06, 1.4, '(b)',
         ha='center',va='center',
         transform = axA1.transAxes, color = expcolor, fontsize=abcd_font_size)

axA2= fig.add_subplot(grid[2:4, 3:6])
#axA2.axvline(x = 5 , ymin=-0.2, ymax = 1.2, ls = '--', 
#            linewidth=trace_width, color=mycolors[6])
axA2.plot(t_plot_Goldbeter_osc, b_trace_osc,linewidth=trace_width, color=goldcolor,
         label='Martiel 1987')
ax3 = axA2.twinx()
ax3.plot(t_plot_Maeda_osc, cAMPi_trace_osc,linewidth=trace_width,color=maedacolor,
         label='Maeda 2004')
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_ylabel('Coupled direct and indirect\n'+r'\n negative feedback model, $cAMP_{i}$, A.U.', color=maedacolor,
                fontsize=label_font_size-3)
ax3.set_ylim([-55,260])

axA2.plot(t_plot_Gregor_osc,gregor_campCyto_trace_osc,linewidth=trace_width, 
         color = gregorcolor,label='Gregor 2010')
axA2.plot(t_plot_Sgro_osc, A_trace_plot_osc,linewidth=trace_width,
         color= sgrocolor, label='Sgro 2015')
axA2.plot(t_plot_Kamino_osc, y_trace_osc, linewidth=trace_width, 
         ls = '--',color=kaminocolor, label='Kamino 2017')
axA2.set_title('Simulation',color = 'k', fontsize = title_font_size)
axA2.text(-0.06, 1.15, '(d)',
         ha='center',va='center',transform = axA2.transAxes, 
         color = simcolor, fontsize=abcd_font_size)

axA2.set_ylim([-0.4,1.85])# oscillations
axA2.set_xlim([0,6])
axA2.tick_params(axis='both', which='major', labelsize=tick_font_size)

axA01.text(-0.13,0.5,' FRET, A.U.', ha='center',va='center',rotation=90,
          transform = axA01.transAxes, color = 'k', fontsize=label_font_size)
axA02.text(-0.13,0.5,r'$cAMP_{i}$, A.U.', ha='center',va='center',rotation=90,
          transform = axA02.transAxes, color = 'k', fontsize=label_font_size)
axA02.text(1.05,-0.27,'Time, A.U.', ha='center',va='center',
          transform = axA02.transAxes, color = 'k', fontsize=label_font_size)
grid.tight_layout(fig,rect=[0, 0, 1, 1],pad = 10)
fig.tight_layout()

# Fig 3 single cell step ramp
# experimental data
my_dir = r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/'
Sgro2015Figure3excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure3')
# load simulation data
import_npz('../model_outputs/single_cell_StepRamp_200320.npz',globals())

# axB0 = fig.add_subplot(grid[0, 0])
# axB0 = fig.add_axes([0.124, 0.72, 0.215, 0.1])
axB0 = fig.add_subplot(grid[6, 0:2])
axB0.plot(Sgro2015Figure3excel["Ramp Input (min Time)"],Sgro2015Figure3excel["Ramp Input (nM cAMP)"],
                              color='k', linewidth = trace_width)
axB0.set_title('Experiment',color='k', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
axB0.set_ylabel(r'$cAMP_{e}$'+'\n(nM)',fontsize=sublabel_font_size)
#axB0.text(-16, 0.5, r'$cAMP_{e}$(nM)', ha='center',va='center',rotation ='vertical',
#         color = 'k', fontsize=label_font_size)
# axB0.set_xlabel('Time, A.U.',fontsize=label_font_size)
axB0.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB0.set_xlim([0,80])
axB0.set_ylim([-0.1,1.1])
axB0.xaxis.set_major_formatter(plt.NullFormatter()) # hide x axis
axB0.text(-0.06, 1.5, '(a)', ha='center',va='center',
         transform = axB0.transAxes, color = expcolor, fontsize=abcd_font_size)

axB1= fig.add_subplot(grid[7:9, 0:2])
axB1.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 1 FRET Trace"],
                               color='k', linewidth = trace_width)
axB1.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 2 FRET Trace"],
                               color='grey', linewidth = trace_width)
axB1.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 3 FRET Trace"],
                               color='lightgrey', linewidth = trace_width)
axB1.axvspan(10, 30, alpha=0.2, color=expcolor); axB1.axvspan(50, 70, alpha=0.2, color='b')
axB1.set_ylabel('FRET, A.U.',fontsize=sublabel_font_size)
axB1.set_xlabel('Time (min)',fontsize=sublabel_font_size)
axB1.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB1.set_ylim([-0.1,0.6]); axB1.set_xlim([0,80])

axB2= fig.add_subplot(grid[6:8,2:4])
axB2.plot(t_plot_Goldbeter, b_trace, color=goldcolor,linewidth=trace_width)
axB2.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB2.set_title('Receptor\n desensitization',color = goldcolor,  fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
axB2.axvspan(2, 6, alpha=0.3, color=simcolor); axB2.axvspan(10, 14, alpha=0.3, color='g')
axB2.set_ylim([-0.2,1.2]); axB2.set_xlim([0,16])
axB2.text(-0.08, 1.2, '(b)', ha='center',va='center',
         transform = axB2.transAxes, color = simcolor, fontsize=abcd_font_size)

axB3= fig.add_subplot(grid[6:8, 4:6])
axB3.plot(t_plot_Maeda, cAMPi_trace, color= maedacolor,linewidth=trace_width)
axB3.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB3.set_title('Coupled direct and\n indirect feedback', color = maedacolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
axB3.axvspan(2,6, alpha=0.3, color=simcolor); axB3.axvspan(10,14, alpha=0.3, color='g')
axB3.set_ylim([-0.2,1.2]);  axB3.set_xlim([0,16])
# axB3.text(-0.12, 1.06, 'C', ha='center',va='center',
#          transform = axB3.transAxes, color = simcolor, fontsize=abcd_font_size)

axB4= fig.add_subplot(grid[9:, 2:4])
axB4.plot(t_plot_Sgro, A_trace_plot, color=sgrocolor,linewidth=trace_width)
axB4.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB4.set_title('Interlocking positive-\n negative feedback', color = sgrocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
axB4.axvspan(2,6, alpha=0.3, color=simcolor); axB4.axvspan(10,14, alpha=0.3, color='g')
axB4.set_ylim([-0.4,1.2]);  axB4.set_xlim([0,16])
# axB4.text(-0.12, 1.06, 'D', ha='center',va='center',
#          transform = axB4.transAxes, color = simcolor, fontsize=abcd_font_size)

axB5= fig.add_subplot(grid[9:, 4:6])
axB5.plot(t_plot_Kamino, y_trace, color=kaminocolor,linewidth=trace_width)
axB5.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB5.set_title('IFFL', color = kaminocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
axB5.axvspan(2,6, alpha=0.3, color= simcolor); axB5.axvspan(10,14, alpha=0.3, color='g')
axB5.set_ylim([-0.4,1.2]);  axB5.set_xlim([0,16])
# axB5.text(-0.12, 1.06, 'E', ha='center',va='center',
#          transform = axB5.transAxes, color = simcolor, fontsize=abcd_font_size)

axB4.text(-0.09,1.26, r'$cAMP_{i}$', ha='center',va='center',rotation=90,
          transform = axB4.transAxes, color = 'k', fontsize=label_font_size)
axB4.text(1.1,-0.4, 'Time, A.U.', ha='center',va='center',
          transform = axB4.transAxes, color = 'k', fontsize=label_font_size)
axB0.text(-0.29, 1.9, 'B',
         ha='center',va='center',
         transform = axB0.transAxes, color = 'k', fontsize=ABCD_font_size)

####################Fig 5 single cell FCD
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


annotate_size = 12
legend_font_size = 11
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Gold)+1))


axC0u = fig.add_axes([0.635,0.78,0.112,0.085],xticks=[0,10,20,30],yticks=[0,2])
# axC0u = fig.add_subplot(grid[0:1, 7:9], xticks=[0,10,20,30],yticks=[0,2])
axC0u.plot(t_plot_Kamino_B,signal_trace, 'k',linewidth=trace_width)
axC0u.set_xlim([0,30]); axC0u.set_ylim([-0.5,3.5]); 
axC0u.tick_params(axis='both', which='major', labelsize=tick_font_size)
axC0u.text(-0.32 , 1.25, '(a)', ha='center',va='center',
     transform = axC0u.transAxes, color = expcolor, fontsize=abcd_font_size)

axC0u.text(-0.28, 0.5, r'$cAMP_{e}$'+'\n input',ha='center',va='center', rotation='vertical',
     transform = axC0u.transAxes, color = 'k', fontsize=sublabel_font_size)

axC0u.annotate(s='', xy=(11,-0.5), xytext=(11,1.2), arrowprops={'arrowstyle':'<->', 'lw':2.5})
axC0u.text(0.5, 0.65, 'Prime\nconc.',ha='center',va='center', 
     transform = axC0u.transAxes, color = 'k', fontsize=annotate_size)

axC0u.annotate(s='', xy=(21,0.7), xytext=(21,3.2), arrowprops={'arrowstyle':'<->', 'lw':2.5})
axC0u.text(0.86, 0.46, 'Fold\nchange',ha='center',va='center', 
     transform = axC0u.transAxes, color = 'k', fontsize=annotate_size)

axC0l = fig.add_axes([0.635,0.64,0.112,0.085],xticks=[0,10,20,30],yticks=[0,0.2,0.4])
# axC0l = fig.add_subplot(grid[2:4,7:9],xticks=[0,10,20,30],yticks=[0,0.2,0.4])
axC0l.plot(t_plot_Kamino_B, y_trace_B, linewidth=trace_width)
axC0l.set_xlim([0,30]); axC0l.set_ylim([-0.05,0.4]); 
axC0l.tick_params(axis='both', which='major', labelsize=tick_font_size)
axC0l.set_xlabel('Time, A.U.', size=sublabel_font_size)
axC0l.text(-0.28, 0.5, r'$cAMP_{i}$'+'\n response',ha='center',va='center', rotation='vertical',
      transform = axC0l.transAxes, color = 'k', fontsize=sublabel_font_size)

axC0l.annotate('', xy=(23,0.02), xytext=(23,0.18),rotation = 'vertical',
            arrowprops={'arrowstyle': '<->','lw':2.5}, va='center')
axC0l.text(0.75, 0.72, 'Second peak\nprominence',ha='center',va='center', 
      transform = axC0l.transAxes, color = 'k', fontsize=annotate_size)

axC0 = fig.add_axes([0.79,0.64,0.115,0.22],yticks=[1.5,2,2.5])
# axC0 = fig.add_subplot(grid[0:4,9:],yticks=[1.5,2,2.5])
#for i in range(len(z0First_space_Sgro)):
#    ax3.plot(FC_space_Sgro_noise,PkPrm_Sgro_mean_noise[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Sgro_noise[i]))
#plt.errorbar(Kamino_FCD["FC_100pM"], Kamino_FCD["100pM mean"], yerr=None, xerr=None, fmt='', ecolor=None, elinewidth=None, capsize=None, barsabove=False, lolims=False, uplims=False, xlolims=False, xuplims=False, errorevery=1, capthick=None, *, data=None, **kwargs)
axC0.errorbar(Kamino_FCD["FC_100pM"], Kamino_FCD["100pM mean"], yerr=Kamino_FCD["100pM SD"], xerr=None, color='SteelBlue', linewidth=trace_width, label='100pM', ecolor='SteelBlue', elinewidth=trace_width,capsize=5,capthick=2)
axC0.errorbar(Kamino_FCD["FC_1nM"], Kamino_FCD["1nM mean"], yerr=Kamino_FCD["1nM SD"], xerr=None, color='SkyBlue', linewidth=trace_width, label='1nM', ecolor='SkyBlue', elinewidth=trace_width,capsize=5,capthick=2)
axC0.errorbar(Kamino_FCD["FC_3nM"], Kamino_FCD["3nM mean"], yerr=Kamino_FCD["3nM SD"], xerr=None, color='turquoise', linewidth=trace_width, label='3nM', ecolor='turquoise', elinewidth=trace_width,capsize=5,capthick=2)
axC0.errorbar(Kamino_FCD["FC_10nM"], Kamino_FCD["10nM mean"], yerr=Kamino_FCD["10nM SD"], xerr=None, color='cyan', linewidth=trace_width, label='10nM', ecolor='cyan', elinewidth=trace_width,capsize=5,capthick=2)
axC0.set_ylim([1.5,2.8])
# axC0.set_ylabel( 'Response Amplitude, A.U.',fontsize=tick_font_size)
axC0.set_xlabel(r'$cAMP_{e}$'+' fold change',fontsize=sublabel_font_size)
axC0.set_xscale('log')
axC0.tick_params(axis='both', which='major', labelsize=tick_font_size)
axC0.set_title('Experiment', color = 'k',fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = axC0.legend();
axC0.legend( frameon=False,loc='bottom center',ncol=1,prop={'size': legend_font_size})
#axC0.text(-0.1, 1.15, 'B', ha='center',va='center',
#         transform = axC0.transAxes, color = 'b', fontsize=abcd_font_size)
axC0.text(-0.23, 0.53, 'Second peak prominence, A.U.', ha='center',va='center', rotation ='vertical',
          transform = axC0.transAxes, fontsize=label_font_size-2)
#
axC1= fig.add_subplot(grid[5:7, 7:9])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Gold)+1) )
for i in range(len(z0First_space_Gold)):
    axC1.plot(FC_space_Gold,PkPrm_Gold[i,:], color=colors[i],
             linewidth=trace_width,label='Priming Conc. '+str(int(z0First_space_Gold[i]*10))+' , A.U.')
# axC1.set_ylim([-0.1,1.2]); 
axC1.set_xlim([FC_space_Gold[0],FC_space_Gold[-1]])
#axC1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#axC1.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
axC1.set_xscale('log')
axC1.tick_params(axis='both', which='major', labelsize=tick_font_size)
# axC1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axC1.set_title('Receptor\n desensitization', color= goldcolor,
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
axC1.text(-0.32, 1.15, '(b)', ha='center',va='center',
         transform = axC1.transAxes, color = simcolor, fontsize=abcd_font_size)
#leg = axC1.legend();
#axC1.legend( frameon=False,loc='bottom center',ncol=1,prop={'size': legend_font_size})

axC2= fig.add_subplot(grid[5:7,9:])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Loomis)+1) )
for i in range(len(z0First_space_Loomis)):
    axC2.plot(FC_space_Loomis,PkPrm_Loomis[i,:],color=colors[i], linewidth=trace_width,label='Prime Conc.'+str(z0First_space_Loomis[i]))

# axC2.set_ylim([0,0.6]); 
axC2.set_xlim([FC_space_Loomis[0],FC_space_Loomis[-1]])
axC2.set_xscale('log')
axC2.tick_params(axis='both', which='major', labelsize=tick_font_size)
axC2.set_title('Coupled direct and\n indirect feedback', color = maedacolor,
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# axC2.text(-0.1, 1.18, 'C', ha='center',va='center',
#          transform = axC2.transAxes, color = simcolor, fontsize=abcd_font_size)

# Sgro with noise
axC5= fig.add_subplot(grid[8:10, 7:9])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Sgro)+1) )
for i in range(len(z0First_space_Sgro_noise)):
    axC5.plot(FC_space_Sgro_noise, PkPrm_Sgro_mean_noise[i,:],'o-', color = colors[i], lw = trace_width, 
             ms = 4,  label='Prime Conc.'+str(z0First_space_Sgro_noise[i]))
    axC5.errorbar(FC_space_Sgro_noise, PkPrm_Sgro_mean_noise[i,:], yerr=PkPrm_Sgro_se_noise[i,:],
                 fmt = 'o', color=colors[i], ecolor= colors[i], elinewidth=trace_width, capsize=5, capthick=2)
axC5.set_ylim([-0.1,1]);
axC5.set_xlim([FC_space_Sgro[0],FC_space_Sgro[-1]])
axC5.set_xscale('log')
axC5.tick_params(axis='both', which='major', labelsize=tick_font_size)
axC5.set_title('Interlocking positive-\n negative feedback', color=sgrocolor,
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# axC5.text(-0.1, 1.18, 'D', ha='center',va='center',
#          transform = axC5.transAxes, color = simcolor, fontsize=abcd_font_size)

axC4= fig.add_subplot(grid[8:10, 9:])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Kamino)+1) )
for i in range(len(z0First_space_Kamino)):
    axC4.plot(FC_space_Kamino,PkPrm_Kamino[i,:],color=colors[i], linewidth=trace_width,label='Prime Conc.'+str(z0First_space_Kamino[i]))
axC4.set_ylim([0,1.1]); 
axC4.set_xlim([FC_space_Kamino[0],FC_space_Kamino[-1]])    
axC4.set_xscale('log')
axC4.tick_params(axis='both', which='major', labelsize=tick_font_size)
axC4.set_title('IFFL',  color=kaminocolor,
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# axC4.text(-0.1, 1.18, 'E', ha='center',va='center',
#          transform = axC4.transAxes, color = simcolor, fontsize=abcd_font_size)
#leg = axC4.legend();
#axC4.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})
axC4.legend(bbox_to_anchor=(0.75,-0.4),ncol=2,prop={'size': 12})

axC5.text(-0.26,1.26,'Second peak prominence', ha='center',va='center',rotation=90,
          transform = axC5.transAxes, color = 'k', fontsize=label_font_size)
axC5.text(1.12,-0.32,  r'$cAMP_{e}$'+' fold change', ha='center',va='center',
          transform = axC5.transAxes, color = 'k', fontsize=label_font_size)

axC0u.text(-0.32, 1.9, 'C',
         ha='center',va='center',
         transform = axC0u.transAxes, color = 'k', fontsize=ABCD_font_size)



grid.tight_layout(fig,rect=[0, 0, 1, 1],pad = 2)
fig.tight_layout()
