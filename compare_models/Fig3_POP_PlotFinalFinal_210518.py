# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 14:53:03 2020

@author: cqhuyan
plotting script for figure 3- pop behaviors

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

expcolor = 'black' #(colors[5,:]+(255-colors[5,:])*0)/255
simcolor = 'gray'# (colors[6,:]+(255-colors[6,:])*0)/255

#%% plotting font sizes
ABCD_font_size = 25
abcd_font_size = 25
label_font_size=16
title_font_size = 16
sublabel_font_size = 14
trace_width=2
tick_font_size=14
text_size = 10.3

fig = plt.figure(figsize=(21,11))
grid = plt.GridSpec(13, 14, wspace=0.8, hspace=0.9)

# population add cAMP, with single cell noise
# Experimental data
my_dir = r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/'
Sgro2015Figure6excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',
                                     sheet_name='Figure6')
# load saved npz output file 
import_npz('../model_outputs/Fig8_pop_add_cAMP_042320.npz',globals())

# experiment
axA01= fig.add_subplot(grid[0:2, 0:3])
axA01.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["Low External cAMP Mean Trace"],
                                color = 'k', linewidth=trace_width)
axA01.axvspan(60, 120, alpha=0.2, color= 'grey')
axA01.set_ylim([-0.1,0.6]); axA01.set_xlim([0,120])
axA01.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA01.text(0.71,0.82,' Low External cAMP, \n 5-10nM', horizontalalignment='center',verticalalignment='center',
          transform = axA01.transAxes, color = 'k', fontsize= text_size)
axA01.set_title('Experiment',size=title_font_size)

axA02= fig.add_subplot(grid[2:4, 0:3])
axA02.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["Intermediate External cAMP Mean Trace"],
                                color = 'k', linewidth=trace_width)
axA02.axvspan(60, 120, alpha=0.2, color = 'grey')
axA02.set_ylim([-0.1,0.6]); axA02.set_xlim([0,120])
axA02.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
axA02.text(0.71,0.82,' Intermediate External \n cAMP, 10-20nM', horizontalalignment='center',verticalalignment='center',
          transform = axA02.transAxes, color = 'k', fontsize=text_size)
axA02.set_ylabel('FRET Signal, A.U.', size=sublabel_font_size)

axA03= fig.add_subplot(grid[4:6, 0:3])
axA03.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["High External cAMP Mean Trace"],
                                color = 'k', linewidth=trace_width)
axA03.axvspan(60, 120, alpha=0.2, color='grey')
axA03.set_ylim([-0.1,0.6]);  axA03.set_xlim([0,120])
axA03.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA03.text(0.71,0.75,' High External cAMP, \n 100nM', horizontalalignment='center',verticalalignment='center',
          transform = axA03.transAxes, color = 'k', fontsize=text_size)
# axA03.set_xlabel('Time (min)', size=sublabel_font_size)
# axA01.text(-0.06 , 1.2, '(a)',
#          horizontalalignment='center',verticalalignment='center',
#          transform = axA01.transAxes, color = expcolor, fontsize=abcd_font_size)
# Goldbeter 1987
axA11= fig.add_subplot(grid[0:2, 3:6])
for count in range(5):
    axA11.plot(t_plot_Goldbeter,b_traces_norm[0,count,:],
              color='darkgrey',alpha=0.5, linewidth=trace_width-1)
axA11.plot(t_plot_Goldbeter,b_traces_norm_norm_mean[0,:], color= goldcolor,linewidth=trace_width)
axA11.text(0.73,0.88,r'Low $cAMP_{e}$ Input', ha='center',va='center',
      transform = axA11.transAxes, color = 'k', fontsize= text_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ Input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
axA11.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA11.axvspan(15, 30, alpha=0.2, color= simcolor)
axA11.set_xlim([0,30]); axA11.set_ylim([-0.35,1.75])

# axA11.text(-0.06 , 1.2, '(b)',
#           horizontalalignment='center',verticalalignment='center',
#           transform = axA11.transAxes, color = simcolor, fontsize=abcd_font_size)
axA11.set_title('Receptor\n Desensitization',color = goldcolor, size=title_font_size)

axA12= fig.add_subplot(grid[2:4, 3:6])
for count in range(5):
    axA12.plot(t_plot_Goldbeter, b_traces_norm[1,count,:],
              color='darkgrey',alpha=0.5, linewidth=trace_width-1)
axA12.plot(t_plot_Goldbeter, b_traces_norm_norm_mean[1,:], color=goldcolor,linewidth=trace_width)
    

axA12.text(0.68,0.88,r'Intermediate $cAMP_{e}$'+' Input', ha='center',va='center',
      transform = axA12.transAxes, color = 'k', fontsize= text_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ Input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
axA12.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA12.axvspan(15, 30, alpha=0.2, color=simcolor)
axA12.set_xlim([0,30]); axA12.set_ylim([-0.35,1.75])

axA13= fig.add_subplot(grid[4:6, 3:6])
for count in range(5):
    axA13.plot(t_plot_Goldbeter, b_traces_norm[2,count,:],
              color='darkgrey',alpha=0.5, linewidth=trace_width-1)
axA13.plot(t_plot_Goldbeter, b_traces_norm_norm_mean[2,:], color=goldcolor,linewidth=trace_width)
axA13.text(0.73,0.88,r'High $cAMP_{e}$ Input', ha='center',va='center',
      transform = axA13.transAxes, color = 'k', fontsize= text_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ Input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
axA13.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
axA13.axvspan(15, 30, alpha=0.2, color=simcolor)
axA13.set_xlim([0,30]); axA13.set_ylim([-0.35,1.75])

# Maeda 2004
axA21= fig.add_subplot(grid[0:2, 6:9])
for count in range(5):
    axA21.plot(t_plot_Maeda_short,cAMPi_traces_norm[0,count,:],
              color='darkgrey',alpha=0.5, linewidth=trace_width-1)
axA21.plot(t_plot_Maeda_short, cAMPi_traces_norm_mean[0,:], color=maedacolor,linewidth=trace_width)
axA21.text(0.73,0.88,r'Low $cAMP_{e}$ Input', ha='center',va='center',
      transform = axA21.transAxes, color = 'k', fontsize= text_size)
axA21.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA21.axvspan(15, 30, alpha=0.2, color = simcolor)
axA21.set_xlim([0,30]); axA21.set_ylim([0.1,0.95])

# axA21.text(-0.06 , 1.2, '(c)',
#           horizontalalignment='center',verticalalignment='center',
#           transform = axA21.transAxes, color = simcolor, fontsize=abcd_font_size)
# axA21.set_title('Coupled Direct-indirect negative feedback',color = maedacolor, size=title_font_size)
axA21.set_title('CDINFB',color = maedacolor, size=title_font_size)

axA22= fig.add_subplot(grid[2:4, 6:9])
for count in range(5):
    axA22.plot(t_plot_Maeda_short,cAMPi_traces_norm[1,count,:],
              color='darkgrey',alpha=0.5, linewidth=trace_width-1)
axA22.plot(t_plot_Maeda_short, cAMPi_traces_norm_mean[1,:], 
          color= maedacolor,linewidth=trace_width)
axA22.text(0.68,0.88,r'Intermediate $cAMP_{e}$ Input', ha='center',va='center',
          transform = axA22.transAxes, color = 'k', fontsize= text_size)
axA22.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA22.axvspan(15, 30, alpha=0.2, color=simcolor)
axA22.set_xlim([0,30]); axA22.set_ylim([0.1,0.95])

axA23= fig.add_subplot(grid[4:6, 6:9])
for count in range(5):
    axA23.plot(t_plot_Maeda_short,cAMPi_traces_norm[2,count,:],
              color= 'darkgrey',alpha=0.5, linewidth=trace_width-1)
axA23.plot(t_plot_Maeda_short, cAMPi_traces_norm_mean[2,:], 
          color=maedacolor,linewidth=trace_width)
axA23.text(0.76,0.88,r'High $cAMP_{e}$ Input', ha='center',va='center',
      transform = axA23.transAxes, color = 'k', fontsize= text_size)
axA23.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA23.axvspan(15, 30, alpha=0.2, color = simcolor)
axA23.set_xlim([0,30]); axA23.set_ylim([0.1,0.95])

# Gregor 2010
axA31= fig.add_subplot(grid[7:9, 0:3])
axA31.plot(t_plot_Gregor[1:],campCyto_traces[0,1:] , color=gregorcolor,linewidth=trace_width)  
for count in range(5):
    # campCyto_traces_single_cell[0,count,:] = campCyto_traces_single_cell[0,count,:]/np.amax(campCyto_traces_single_cell[0,count,:])
    axA31.plot(t_plot_Gregor,campCyto_traces_single_cell[0,count,:],
              color='darkgrey',alpha=0.5, linewidth=trace_width-1)
 
axA31.text(0.73,0.88,r'Low $cAMP_{e}$ Input', ha='center',va='center',
      transform = axA31.transAxes, color = 'k', fontsize= text_size)
axA31.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA31.axvspan(15, 30, alpha=0.2, color= simcolor)
axA31.set_xlim([0,30]); axA31.set_ylim([-0.35,1.75])

# axA31.text(-0.06 , 1.2, '(d)', ha='center',va='center',
#          transform = axA31.transAxes, color = simcolor, fontsize=abcd_font_size)
axA31.set_title('Phase Oscillator',color = gregorcolor, size=title_font_size)

axA32= fig.add_subplot(grid[9:11, 0:3])
axA32.plot(t_plot_Gregor[1:],campCyto_traces[1,1:], color=gregorcolor,linewidth=trace_width)
for count in range(5):
    # campCyto_traces_single_cell[1,count,:] = campCyto_traces_single_cell[1,count,:]/np.amax(campCyto_traces_single_cell[1,count,:])
    axA32.plot(t_plot_Gregor,campCyto_traces_single_cell[1,count,:],
              color='darkgrey',alpha=0.5, linewidth=trace_width-1)
axA32.text(0.68,0.88,r'Intermediate $cAMP_{e}$'+' Input', ha='center',va='center',
      transform = axA32.transAxes, color = 'k', fontsize= text_size)
axA32.tick_params(grid_linewidth =tick_font_size, labelsize = tick_font_size)
axA32.axvspan(15, 30, alpha=0.2, color = simcolor)
axA32.set_xlim([0,30]); axA32.set_ylim([-0.35,1.75])

axA33= fig.add_subplot(grid[11:, 0:3])
axA33.plot(t_plot_Gregor[1:],campCyto_traces[2,1:], color= gregorcolor, linewidth=trace_width)
for count in range(5):
    # campCyto_traces_single_cell[2,count,:] = campCyto_traces_single_cell[2,count,:]/np.amax(campCyto_traces_single_cell[2,count,:])
    axA33.plot(t_plot_Gregor,campCyto_traces_single_cell[2,count,:],
              color='darkgrey',alpha=0.5, linewidth=trace_width-1)
axA33.text(0.73,0.88,r'High $cAMP_{e}$ Input', ha='center',va='center',
      transform = axA33.transAxes, color = 'k', fontsize= text_size)
axA33.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA33.axvspan(15, 30, alpha=0.2, color=simcolor)
axA33.set_xlim([0,30]); axA33.set_ylim([-0.35,1.75])

# Sgro 2015
axA41= fig.add_subplot(grid[7:9, 3:6])
for count in range(5):
    axA41.plot(t_plot_Sgro,A_traces_single_cell[0,count,:],
              color= 'darkgrey' ,alpha=0.3, linewidth=trace_width-1)
axA41.plot(t_plot_Sgro,A_traces[0,:], color= sgrocolor,linewidth=trace_width)
axA41.text(0.73,0.88,r'Low $cAMP_{e}$ Input', ha='center',va='center',
      transform = axA41.transAxes, color = 'k', fontsize= text_size)
axA41.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA41.axvspan(15, 30, alpha=0.2, color=simcolor)
axA41.set_xlim([0,30]); axA41.set_ylim([-0.35,1.75])

# axA41.text(-0.06 , 1.2, '(e)',ha='center',va='center',
#           transform = axA41.transAxes, color = simcolor, fontsize=abcd_font_size)
# axA41.set_title('Interlocking positive-\n negative feedback',color = sgrocolor, size=title_font_size)
axA41.set_title('IPNFB',color = sgrocolor, size=title_font_size)

axA42= fig.add_subplot(grid[9:11, 3:6])
for count in range(5):
    axA42.plot(t_plot_Sgro,A_traces_single_cell[1,count,:],
              color='darkgrey',alpha=0.3, linewidth=trace_width-1)
axA42.plot(t_plot_Sgro,A_traces[1,:], color=sgrocolor,linewidth=trace_width)
    

axA42.text(0.68,0.88,r'Intermediate $cAMP_{e}$'+' Input', ha='center',va='center',
      transform = axA42.transAxes, color = 'k', fontsize= text_size)
axA42.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA42.axvspan(15, 30, alpha=0.2, color=simcolor)
axA42.set_xlim([0,30]); axA42.set_ylim([-0.35,1.75])

axA43= fig.add_subplot(grid[11:, 3:6])
for count in range(5):
    axA43.plot(t_plot_Sgro,A_traces_single_cell[2,count,:],
              color='darkgrey',alpha=0.3, linewidth=trace_width-1)
axA43.plot(t_plot_Sgro,A_traces[2,:], color=sgrocolor,linewidth=trace_width)
axA43.text(0.73,0.88,r'High $cAMP_{e}$ Input', ha='center',va='center',
      transform = axA43.transAxes, color = 'k', fontsize= text_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ Input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
axA43.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA43.axvspan(15, 30, alpha=0.2, color=simcolor)
axA43.set_xlim([0,30]); axA43.set_ylim([-0.35,1.75])

# Kamino 2017
axA51= fig.add_subplot(grid[7:9, 6:9])
for count in range(5):
    axA51.plot(t_plot_Kamino,y_traces_norm[0,count,:],
              color='darkgrey',alpha=0.5, linewidth=trace_width-1)
axA51.plot(t_plot_Kamino,y_traces_norm_mean[0,:], color= kaminocolor,linewidth=trace_width)
axA51.text(0.73,0.88,r'Low $cAMP_{e}$ Input', ha='center',va='center',
      transform = axA51.transAxes, color = 'k', fontsize= text_size)
axA51.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA51.axvspan(15, 30, alpha=0.2, color=simcolor)
axA51.set_xlim([0,30]); axA51.set_ylim([-0.35,1.75])

# axA51.text(-0.06 , 1.2, '(f)', ha='center',va='center',
#           transform = axA51.transAxes, color = simcolor, fontsize=abcd_font_size)
axA51.set_title('IFFL',color = kaminocolor, size=title_font_size)

axA52= fig.add_subplot(grid[9:11, 6:9])
for count in range(5):
    axA52.plot(t_plot_Kamino,y_traces_norm[1,count,:],
              color='darkgrey',alpha=0.5, linewidth=trace_width-1)
axA52.plot(t_plot_Kamino,y_traces_norm_mean[1,:], color= kaminocolor, linewidth=trace_width)
axA52.text(0.73,0.75,r'Intermediate $cAMP_{e}$'+'\n Input', ha='center',va='center',
      transform = axA52.transAxes, color = 'k', fontsize= text_size)
axA52.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA52.axvspan(15, 30, alpha=0.2, color=simcolor)
axA52.set_xlim([0,30]); axA52.set_ylim([-0.35,1.75])

axA53= fig.add_subplot(grid[11:, 6:9])
for count in range(5):
    axA53.plot(t_plot_Kamino,y_traces_norm[2,count,:],
              color='darkgrey',alpha=0.5, linewidth=trace_width-1)
axA53.plot(t_plot_Kamino,y_traces_norm_mean[2,:], color= kaminocolor, linewidth=trace_width)
axA53.text(0.73,0.88,r'High $cAMP_{e}$ Input', ha='center',va='center',
     transform = axA53.transAxes, color = 'k', fontsize= text_size)
axA53.tick_params(grid_linewidth = tick_font_size, labelsize = tick_font_size)
axA53.axvspan(15, 30, alpha=0.2, color=simcolor)
axA53.set_xlim([0,30]); axA53.set_ylim([-0.35,1.75]) 

axA33.text(-0.15,1.65, r'$cAMP_{i}$', ha='center',va='center',rotation=90,
          transform = axA33.transAxes, color = 'k', fontsize=label_font_size)
axA33.text(1.75,-0.5, 'Time, A.U.', ha='center',va='center',
          transform = axA33.transAxes, color = 'k', fontsize=label_font_size)


# pop firing rate

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

# Gregor 2010 from parallel computing outputs
Gregor_OUT_noise = np.load(OUT_path +'pop_fire_rate_Gregor_OUT_191026_noise0.002.npz')
k_arr_Gregor_noise =  Gregor_OUT_noise['k_arr']
rho_arr_Gregor_noise =  Gregor_OUT_noise['rho_arr']
pop_rate_Gregor_noise = Gregor_OUT_noise['pop_rate_Gregor']

Gregor_OUT = np.load(OUT_path +'pop_fire_rate_Gregor_OUT_191027_noise0.npz')
k_arr_Gregor =  Gregor_OUT['k_arr']
rho_arr_Gregor =  Gregor_OUT['rho_arr']
pop_rate_Gregor = Gregor_OUT['pop_rate_Gregor']

# Sgro 2015
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

# Kamino 2017 from parallel computing outputs
Kamino_OUT = np.load(OUT_path +'pop_FR_Kamino_200401_dt0.001_noise_0_PrmLen_25PkFindThr0.03.npz')
rho_arr_Kamino =  Kamino_OUT['rho_arr']
gamma_arr_Kamino =  Kamino_OUT['gamma_arr']
pop_rate_Kamino = Kamino_OUT['pop_rate_Kamino']

Kamino_OUT_noise = np.load(OUT_path +'pop_FR_Kamino_200401_dt0.001_noise_0.01_PrmLen_25PkFindThr0.03.npz')
rho_arr_Kamino_noise =  Kamino_OUT_noise['rho_arr']
gamma_arr_Kamino_noise =  Kamino_OUT_noise['gamma_arr']
pop_rate_Kamino_noise = Kamino_OUT_noise['pop_rate_Kamino']


# fig = plt.figure(figsize=(24,9))
# grid = plt.GridSpec(13, 14, wspace=1.5, hspace=0.8)

# axB0= fig.add_subplot(grid[0:3,10:12])
# axB0 = fig.add_axes([0.67, 0.71, 0.09,0.18])
axB0 = fig.add_axes([0.675, 0.675, 0.11,0.24])
axB0.set_xticks([0,2,4,6,8]); 
axB0.set_xticklabels([1,3,5,7,9],fontsize=tick_font_size)

# axB0.set_yticks([0,1,2,3,4,5,6,7]); 
# axB0.set_yticklabels(['1/2','1/4','1/8','1/16','1/32','1/64','1/128'],fontsize=tick_font_size-3)
axB0.set_yticks([0,2,4,6]); 
axB0.set_yticklabels(['1/2','1/8','1/32','1/128'],fontsize=tick_font_size)

axB0.set_title('Experiment', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
axB0.set_xlabel('Flow Rate (mL/min)', size=sublabel_font_size)
axB0.set_ylabel('Cell Density(mML)', size=sublabel_font_size)
heatmap = axB0.imshow(PopRateExp, cmap='jet') # cmap='jet'
x=[3.5,4.5,5.5,7.5,9.5]
[axB0.axvline(_x, color='white',linewidth=trace_width) for _x in x]
# heatmap.set_clim(0,0.16)
cbar=fig.colorbar(heatmap, ax=axB0,ticks=[0,0.05,0.1,0.15]);
cbar.ax.tick_params(labelsize = tick_font_size) 
cbar.set_label( 'cAMP pulses/min',size=tick_font_size)
#axB0.tick_params(axis='both', which='major', labelsize=tick_font_size)
# axB0.text(-0.33 , 1.55, '(a)',
#          ha='center',va='center',
#          transform = axB0.transAxes, color = expcolor, fontsize=abcd_font_size)


# axB1lower= fig.add_subplot(grid[0:3,12:],xticks=[0,25,50,75,100])
# axB1lower = fig.add_axes([0.82, 0.68, 0.11,0.17],xticks=[0,25,50,75,100])
axB1lower = fig.add_axes([0.84, 0.675, 0.135,0.205],xticks=[0,25,50,75,100])
heatmap = axB1lower.pcolor(kc_arr_Gold, oneoverh_arr_Gold, pop_rate_Gold_noise10.transpose(), cmap='jet') # cmap='jet'
# axB1.set_xscale('log');
axB1lower.set_yscale('log')
heatmap.set_clim(0,1.5)
cbar=fig.colorbar(heatmap, ax=axB1lower, ticks=[0,0.5,1,1.5])
cbar.ax.tick_params(labelsize = tick_font_size) 
axB1lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB1lower.set_title('Receptor\nDesensitization', color = goldcolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# axB1lower.text(-0.25 , 1.2, '(b)',
#           ha='center',va='center',
#           transform = axB1lower.transAxes, color = simcolor, fontsize=abcd_font_size)


# axB2lower= fig.add_subplot(grid[4:7,10:12],xticks=[0,25,50,75,100])
# axB2lower = fig.add_axes([0.672, 0.44, 0.11,0.17], xticks=[0,25,50,75,100])
axB2lower = fig.add_axes([0.675, 0.39, 0.135,0.205],xticks=[0,25,50,75,100])
heatmap = axB2lower.pcolor(gamma_arr_Maeda, rho_arr_Maeda, pop_rate_Maeda_noise1.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.65)
axB2lower.set_yscale('log')
cbar=fig.colorbar(heatmap, ax=axB2lower);cbar.ax.tick_params(labelsize = tick_font_size) 
axB2lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
# axB2lower.set_title('Coupled direct and\nindirect Negative feedback', color= maedacolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
axB2lower.set_title('CDINFB', color= maedacolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# axB2lower.text(-0.25 , 1.2, '(c)',
#           ha='center',va='center',
#           transform = axB2lower.transAxes, color = simcolor, fontsize=abcd_font_size)

# axB3= fig.add_subplot(grid[4:7,12:], xticks=[0,25,50,75,100])
# axB3 = fig.add_axes([0.82, 0.44, 0.11,0.17],xticks=[0,25,50,75,100])
axB3 = fig.add_axes([0.84, 0.39, 0.135,0.205],xticks=[0,25,50,75,100])

heatmap = axB3.pcolor(k_arr_Gregor_noise, rho_arr_Gregor_noise, pop_rate_Gregor_noise.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,1.2)
axB3.set_yscale('log')
cbar=fig.colorbar(heatmap, ax=axB3);cbar.ax.tick_params(labelsize = tick_font_size) 
axB3.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB3.set_title('Phase Oscillator',color= gregorcolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# axB3.text(-0.25 , 1.2, '(d)',
#          ha='center',va='center',
#          transform = axB3.transAxes, color = simcolor, fontsize=abcd_font_size)

# Sgro regular noise (sig = 0.15)
# axB4= fig.add_subplot(grid[8:11,10:12],xticks=[0,0.5,1])
# axB4 = fig.add_axes([0.672, 0.17, 0.11,0.17], xticks=[0,0.5,1])
axB4 = fig.add_axes([0.675, 0.11, 0.135,0.205], xticks=[0,0.5,1])

heatmap = axB4.pcolor(j_arr_Sgro_regular_noise, rho_arr_Sgro_regular_noise, 
                     pop_rate_Sgro_regular_noise.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.65)
axB4.set_yscale('log'); # axB4.set_ylim([10**(-5),10**(-3)]); 
cbar=fig.colorbar(heatmap, ax=axB4);cbar.ax.tick_params(labelsize = tick_font_size) 
axB4.tick_params(axis='both', which='major', labelsize=tick_font_size)
# axB4.set_title('Interlocking positive-\nnegative feedback', color= sgrocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
axB4.set_title('IPNFB', color= sgrocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# axB4.text(-0.25 , 1.2, '(e)',
#          ha='center',va='center',
#          transform = axB4.transAxes, color = simcolor, fontsize=abcd_font_size)
# axB4.set_xticklabels([0,0.25,0.5,0.75,1], rotation=45,fontsize=tick_font_size)


# axB5lower= fig.add_subplot(grid[8:11,12:], xticks=[0,25,50,75,100])
# axB5lower = fig.add_axes([0.82, 0.17, 0.11,0.17],xticks=[0,25,50,75,100])
axB5lower = fig.add_axes([0.84, 0.11, 0.135,0.205],xticks=[0,25,50,75,100])

heatmap = axB5lower.pcolor(gamma_arr_Kamino, rho_arr_Kamino, 
                          pop_rate_Kamino_noise.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.65)
axB5lower.set_yscale('log')
cbar=fig.colorbar(heatmap, ax=axB5lower);cbar.ax.tick_params(labelsize = tick_font_size) 
axB5lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB5lower.set_title('IFFL', color = kaminocolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

axB5lowerin= fig.add_axes([0.895,0.12,0.05,0.086])
heatmap = axB5lowerin.pcolor(gamma_arr_Kamino, rho_arr_Kamino,
                       pop_rate_Kamino_noise.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.65)
axB5lowerin.set_yscale('log'); axB5lowerin.set_xscale('log');
axB5lowerin.set_xticks([]) ; axB5lowerin.set_yticks([]) 
axB5lowerin.spines['bottom'].set_color('white');axB5lowerin.spines['top'].set_color('white')
axB5lowerin.spines['left'].set_color('white');axB5lowerin.spines['right'].set_color('white')

# axB5lower.text(-0.25 , 1.2, '(f)',
#           ha='center',va='center',
#           transform = axB5lower.transAxes, color = simcolor, fontsize=abcd_font_size)

# axA01.text(-0.2, 1.32, 'A',
#          ha='center',va='center',
#          transform = axA01.transAxes, color = 'k', fontsize=ABCD_font_size)
# axB0.text(-0.47, 1.49, 'B',
#          ha='center',va='center',
#          transform = axB0.transAxes, color = 'k', fontsize=ABCD_font_size)
axB4.text(1.36 , -0.22, 'Dilution Rate, A.U.',
         ha='center',va='center',
         transform = axB4.transAxes, color = 'k', fontsize= label_font_size)
axB4.text(-0.32 , 1.2, 'Population Density, A.U.',
         ha='center',va='center',rotation = 90,
         transform = axB4.transAxes, color = 'k', fontsize= label_font_size)

fig.text(0.095, 0.91, 'A',fontsize=ABCD_font_size, ha='center')
fig.text(0.633, 0.91, 'B',fontsize=ABCD_font_size, ha='center')
# fig.text(0.51, 0.045, 'Dilution Rate, A.U.',fontsize=label_font_size, ha='center')
# fig.text(0.08, 0.35, 'Population Density, A.U.',fontsize=label_font_size, va='center', rotation='vertical')

grid.tight_layout(fig,rect=[0, 0, 1, 1],pad = 2)
fig.tight_layout()