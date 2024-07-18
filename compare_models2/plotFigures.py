# plot figures for Chuqiao's paper 

import os, pickle 
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']     

abcd_font_size = 28
label_font_size=16
title_font_size = 18
sublabel_font_size = 22
trace_width=3
tick_font_size=16
letterLabelSize = 32

#%% Figure 2 (Normalizing single spike): pull data 

# experimental data
my_dir = '//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/exp_data/'
Sgro2015Figure1excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure1')

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_singleSpike_data/singleSpike_Gregor.pickle','rb') as f:
    t_Gregor,cAMP_Gregor = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_singleSpike_data/singleSpike_Sgro.pickle','rb') as f:
    t_Sgro,cAMP_Sgro = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_singleSpike_data/singleSpike_Goldbeter.pickle','rb') as f:
    t_Goldbeter,cAMP_Goldbeter = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_singleSpike_data/singleSpike_MaedaLoomis.pickle','rb') as f:
    t_ML,cAMP_ML = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_singleSpike_data/singleSpike_Kamino.pickle','rb') as f:
    t_Kamino,cAMP_Kamino = pickle.load(f)

#%% Figure 2 (Normalizing single spike): plot 

f = plt.figure(figsize=(8,8))
grid = plt.GridSpec(3,2,wspace=0.5,hspace=0.7)

f.text(0.01,0.955,'A',fontsize=letterLabelSize)
ax = f.add_subplot(grid[0,0])
ax.set_title('Experiment:\nAdaptive Spiking',fontsize=title_font_size)
ax.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 1 FRET Trace (1nM)"],
                              linewidth=trace_width,color='k')
ax.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 2 FRET Trace (1nM)"],
                               linewidth=trace_width,color='dimgrey')
ax.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 3 FRET Trace (1nM)"],
                               linewidth=trace_width,color='darkgrey')
ax.text(17,0.5,'1$nm$ cAMP',fontsize=tick_font_size)
ax.vlines(5,-1,1,color='k',linestyle='dashed')
ax.set_xlabel('Time (min)',fontsize=label_font_size)
ax.set_xlim([0,30])
ax.set_ylabel('FRET (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.1,.6])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
f.add_artist(mpatches.FancyBboxPatch((0.033,0.7),0.43,0.268,facecolor='lightgrey',zorder=0,alpha=1,mutation_scale=0.1,linewidth=0))

f.text(0.52,0.955,'B',fontsize=letterLabelSize)
ax = f.add_subplot(grid[0,1])
ax.set_title('Phase Oscillator',fontsize=title_font_size,color=mycolors[2])
ax.plot(t_Gregor,cAMP_Gregor,color='dimgrey',linewidth=2)
ax.arrow(5,1,0,19,color='b',linewidth=2,length_includes_head=True,head_width=1)
ax.arrow(5,20,0,-19,color='b',linewidth=2,length_includes_head=True,head_width=1)
ax.arrow(7,0.5,4.5,0,color='r',linewidth=2,length_includes_head=True,head_width=1,head_length=0.5)
ax.arrow(11.5,0.5,-5,0,color='r',linewidth=2,length_includes_head=True,head_width=1,head_length=0.5)
ax.vlines(6,-2,22,color='k',linestyle='dashed')
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xlim([0,25])
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size)
ax.set_ylim([-2,22])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.add_artist(mpatches.FancyBboxPatch((0.033,0.1),0.935,0.535,edgecolor='k',facecolor='white',zorder=0,mutation_scale=0.1))
f.text(0.01,0.62,'C',fontsize=letterLabelSize)
ax = f.add_subplot(grid[1,0])
ax.set_title('Receptor Desensitization',fontsize=title_font_size,color=mycolors[3])
ax.plot(t_Goldbeter,cAMP_Goldbeter,color='dimgrey',linewidth=2)
ax.vlines(6.941,-100,250,color='k',linestyle='dashed')
ax.arrow(6,10,0,205,color='b',linewidth=2,length_includes_head=True,head_width=1,head_length=15)
ax.arrow(6,205,0,-195,color='b',linewidth=2,length_includes_head=True,head_width=1,head_length=15)
ax.arrow(7.5,25,4.5,0,color='r',linewidth=2,length_includes_head=True,head_width=15,head_length=0.8)
ax.arrow(7.5+4.5,25,-4.5,0,color='r',linewidth=2,length_includes_head=True,head_width=15,head_length=0.8)
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xlim([0,25])
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size,labelpad=-4)
ax.set_ylim([-10,230])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.text(0.52,0.62,'D',fontsize=letterLabelSize)
ax = f.add_subplot(grid[1,1])
ax.set_title('CDINFB',fontsize=title_font_size,color=mycolors[8])
ax.plot(t_ML,cAMP_ML,color='dimgrey',linewidth=2)
ax.vlines(3.83,-100,250,color='k',linestyle='dashed')
ax.arrow(3.1,0.2,0,2.9,color='b',linewidth=2,length_includes_head=True,head_width=0.8,head_length=.25)
ax.arrow(3.1,0.2+2.9,0,-2.9,color='b',linewidth=2,length_includes_head=True,head_width=0.8,head_length=.25)
ax.arrow(5.5,2,2.2,0,color='r',linewidth=2,length_includes_head=True,head_width=0.2,head_length=0.4)
ax.arrow(5.5+2.2,2,-2.2,0,color='r',linewidth=2,length_includes_head=True,head_width=0.2,head_length=0.4)
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xlim([0,20])
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.2,3.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.text(0.01,0.32,'E',fontsize=letterLabelSize)
ax = f.add_subplot(grid[2,0])
ax.set_title('IPNFB',fontsize=title_font_size,color=mycolors[5])
ax.plot(t_Sgro,cAMP_Sgro-np.min(cAMP_Sgro),color='dimgrey',linewidth=2)
ax.vlines(25,-100,250,color='k',linestyle='dashed')
ax.arrow(15,0.8,0,4.2-0.8,color='b',linewidth=2,length_includes_head=True,head_width=3,head_length=.25,zorder=2)
ax.arrow(15,4.2,0,-(4.2-0.8),color='b',linewidth=2,length_includes_head=True,head_width=3,head_length=.25,zorder=2)
ax.arrow(15,0,0,0.7,color='k',linewidth=2,length_includes_head=True,head_width=3,head_length=.25,zorder=2)
ax.arrow(15,0.7,0,-0.7,color='k',linewidth=2,length_includes_head=True,head_width=3,head_length=.25,zorder=2)
ax.arrow(28,0.5,38-28,0,color='r',linewidth=2,length_includes_head=True,head_width=0.3,head_length=2)
ax.arrow(38,0.5,-10,0,color='r',linewidth=2,length_includes_head=True,head_width=0.3,head_length=2)
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xlim([0,80])
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size)
ax.set_ylim([0,4.3])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.text(0.52,0.32,'F',fontsize=letterLabelSize)
ax = f.add_subplot(grid[2,1])
ax.set_title('IFFL',fontsize=title_font_size,color=mycolors[0])
ax.plot(t_Kamino,cAMP_Kamino,linewidth=2,color='dimgrey')
ax.vlines(5.17,-100,250,color='k',linestyle='dashed')
ax.arrow(4,0.07,0,.31-0.07,color='b',linewidth=2,length_includes_head=True,head_width=0.75,head_length=.02,zorder=2)
ax.arrow(4,0.31,0,-(.31-0.07),color='b',linewidth=2,length_includes_head=True,head_width=0.75,head_length=.02,zorder=2)
ax.arrow(4,0,0,.05,color='k',linewidth=2,length_includes_head=True,head_width=0.75,head_length=.02,zorder=2)
ax.arrow(4,0.05,0,-.05,color='k',linewidth=2,length_includes_head=True,head_width=0.75,head_length=.02,zorder=2)
ax.arrow(5.5,0.07,4,0,color='r',linewidth=2,length_includes_head=True,head_width=0.02,head_length=0.5)
ax.arrow(5.5+4,0.07,-4,0,color='r',linewidth=2,length_includes_head=True,head_width=0.02,head_length=0.5)
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xlim([0,20])
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size)
ax.set_ylim([0,0.35])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.add_artist(mpatches.FancyBboxPatch((0.033,0.035),0.935,0.002,edgecolor='k',facecolor='white',zorder=0,mutation_scale=0.1))
f.add_artist(mpatches.FancyArrow(0.12,0.035,0.1,0,color='b',head_width=.015,head_length=0.015,width=0.001))
f.add_artist(mpatches.FancyArrow(0.22,0.035,-0.1,0,color='b',head_width=.015,head_length=0.015,width=0.001))
f.text(0.25,0.028,'Height',fontsize=14)
f.add_artist(mpatches.FancyArrow(0.42,0.035,0.1,0,color='r',head_width=.015,head_length=0.015,width=0.001))
f.add_artist(mpatches.FancyArrow(0.52,0.035,-0.1,0,color='r',head_width=.015,head_length=0.015,width=0.001))
f.text(0.55,0.028,'Width',fontsize=14)
f.add_artist(mpatches.FancyArrow(0.72,0.035,0.1,0,color='k',head_width=.015,head_length=0.015,width=0.001))
f.add_artist(mpatches.FancyArrow(0.82,0.035,-0.1,0,color='k',head_width=.015,head_length=0.015,width=0.001))
f.text(0.85,0.028,'Offset',fontsize=14)
plt.subplots_adjust(top = 0.92, bottom = 0.15, right = 0.98, left = 0.12)
plt.show()

#%% Figure 3 (population add cAMP): pull data

# experimental data
my_dir = '//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/exp_data/'
Sgro2015Figure6excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure6')

# Gregor
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/popAddcAMPscNoise_data/popAddcAMPscNoise_Gregor.pickle','rb') as f:
    ts_Gregor,thetais_Gregor,cAMPes_Gregor,cAMPis_sc_Gregor,cAMPis_Gregor = pickle.load(f)

# Sgro 
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/popAddcAMPscNoise_data/popAddcAMPscNoise_Sgro.pickle','rb') as f:
    ts_Sgro,As_Sgro,As_sc_Sgro = pickle.load(f)

# Goldbeter
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/popAddcAMPscNoise_data/popAddcAMPscNoise_Goldbeter.pickle','rb') as f:
    ts_Gb,ps_Gb,bs_Gb,gs_Gb = pickle.load(f)
    
# Maeda
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/popAddcAMPscNoise_data/popAddcAMPscNoise_Maeda.pickle','rb') as f:
    t_ML,cAMPis_ML = pickle.load(f)

# Kamino
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/popAddcAMPscNoise_data/popAddcAMPscNoise_Kamino.pickle','rb') as f:
    t_Kamino,xs_Kamino,ys_Kamino,zs_Kamino = pickle.load(f)


#%% Figure 3 (Population add cAMP): plot

fig = plt.figure(figsize=(16,12))
grid = plt.GridSpec(7, 3, wspace=0.25, hspace=0.4)

fig.text(0.01,0.96,'A',fontsize=letterLabelSize)
ax= fig.add_subplot(grid[0, 0])
ax.set_title('Experiment',size=title_font_size)
ax.vlines(60,-1,1,color='k',linestyle='dashed')
ax.plot(Sgro2015Figure6excel["Times (min)"],
          Sgro2015Figure6excel["Low External cAMP Mean Trace"],
          color = 'k', linewidth=trace_width)
ax.text(85,0.35,'Low $cAMP_e$ Input\n5-10nM',fontsize=12,horizontalalignment='center')
ax.set_ylim([-0.1,0.6]);
ax.set_xlim([0,120])
ax.set_xticks(np.linspace(0,120,7))
ax.set_xticks(np.linspace(0,120,13),[],minor=1)
ax.set_yticks(np.linspace(0,.5,2))
ax.set_yticks(np.linspace(0,.5,6),[],minor=1)
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax= fig.add_subplot(grid[1,0])
ax.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["Intermediate External cAMP Mean Trace"],
                                color = 'k', linewidth=trace_width)
ax.vlines(60,-1,1,color='k',linestyle='dashed')
ax.text(90,0.35,'Moderate $cAMP_e$\nInput 10-20nM',fontsize=12,horizontalalignment='center')
ax.set_ylim([-0.1,0.6])
ax.set_xlim([0,120])
ax.set_ylabel('FRET (A.U.)',fontsize=label_font_size)
ax.set_xticks(np.linspace(0,120,7))
ax.set_xticks(np.linspace(0,120,13),[],minor=1)
ax.set_yticks(np.linspace(0,.5,2))
ax.set_yticks(np.linspace(0,.5,6),[],minor=1)
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax= fig.add_subplot(grid[2, 0])
ax.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["High External cAMP Mean Trace"],color = 'k', linewidth=trace_width)
ax.vlines(60,-1,1,color='k',linestyle='dashed')
ax.text(90,0.35,'High $cAMP_e$ Input\n100nM',fontsize=12,horizontalalignment='center')
ax.set_ylim([-0.1,0.6])
ax.set_xlim([0,120])
ax.set_xlabel('Time (min)',fontsize=label_font_size)
ax.set_xticks(np.linspace(0,120,7))
ax.set_xticks(np.linspace(0,120,13),[],minor=1)
ax.set_yticks(np.linspace(0,.5,2))
ax.set_yticks(np.linspace(0,.5,6),[],minor=1)
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

fig.add_artist(mpatches.FancyBboxPatch((0.033,0.57),0.28,0.395,facecolor='lightgrey',zorder=0,alpha=1,mutation_scale=0.1,linewidth=0))

fig.text(0.35,0.96,'B',fontsize=letterLabelSize)
ax = fig.add_subplot(grid[0,1])
ax.set_title('Phase Oscillator',color=mycolors[2],fontsize=title_font_size)
for i in range(10):
    ax.plot(ts_Gregor,cAMPis_sc_Gregor[0,i], alpha=0.8, color='grey',linewidth=2)  
ax.plot(ts_Gregor, cAMPis_Gregor[0], alpha=0.8, color=mycolors[2],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'Low $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylim([-0.1,1.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = fig.add_subplot(grid[1,1])
for i in range(10):
    ax.plot(ts_Gregor,cAMPis_sc_Gregor[1,i], alpha=0.8, color='grey',linewidth=2)  
ax.plot(ts_Gregor, cAMPis_Gregor[1], alpha=0.8, color=mycolors[2],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'Moderate $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.1,1.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = fig.add_subplot(grid[2,1])
for i in range(10):
    ax.plot(ts_Gregor,cAMPis_sc_Gregor[2,i], alpha=0.8, color='grey',linewidth=2)  
ax.plot(ts_Gregor, cAMPis_Gregor[2], alpha=0.8, color=mycolors[2],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'High $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylim([-0.1,1.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

fig.text(0.68,0.96,'C',fontsize=letterLabelSize)
ax = fig.add_subplot(grid[0,2])
ax.set_title('IPNFB',color=mycolors[5],fontsize=title_font_size)
for i in range(10):
    ax.plot(ts_Sgro,As_sc_Sgro[0,i], alpha=0.2, color='grey',linewidth=2)  
ax.plot(ts_Sgro, As_Sgro[0], alpha=0.8, color=mycolors[5],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'Low $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylim([-0.3,1.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = fig.add_subplot(grid[1,2])
for i in range(10):
    ax.plot(ts_Sgro,As_sc_Sgro[1,i], alpha=0.2, color='grey',linewidth=2)  
ax.plot(ts_Sgro, As_Sgro[1], alpha=0.8, color=mycolors[5],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'Moderate $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.3,1.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = fig.add_subplot(grid[2,2])
for i in range(10):
    ax.plot(ts_Sgro,As_sc_Sgro[2,i], alpha=0.2, color='grey',linewidth=2)  
ax.plot(ts_Sgro, As_Sgro[2], alpha=0.8, color=mycolors[5],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'High $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylim([-0.3,1.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

fig.add_artist(mpatches.FancyBboxPatch((0.7,0.57),0.265,0.395,facecolor='palegreen',zorder=0,alpha=0.5,mutation_scale=0.1,linewidth=0))

fig.text(0.01,0.44,'D',fontsize=letterLabelSize)
ax = fig.add_subplot(grid[4,0])
ax.set_title('Receptor Desensitization',color=mycolors[7],fontsize=title_font_size)
for i in range(10):
    ax.plot(ts_Gb, bs_Gb[0,i], alpha=0.2, color='grey',linewidth=2)  
ax.plot(ts_Gb, np.mean(bs_Gb[0],axis=0), alpha=0.8, color=mycolors[7],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'Low $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylim([-0.1,1.7])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = fig.add_subplot(grid[5,0])
for i in range(10):
    ax.plot(ts_Gb, bs_Gb[1,i], alpha=0.2, color='grey',linewidth=2)  
ax.plot(ts_Gb, np.mean(bs_Gb[1],axis=0), alpha=0.8, color=mycolors[7],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'Moderate $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.1,1.7])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = fig.add_subplot(grid[6,0])
for i in range(10):
    ax.plot(ts_Gb, bs_Gb[2,i], alpha=0.2, color='grey',linewidth=2)  
ax.plot(ts_Gb, np.mean(bs_Gb[2],axis=0), alpha=0.8, color=mycolors[7],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'High $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylim([-0.1,1.7])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

fig.text(0.35,0.44,'E',fontsize=letterLabelSize)
ax = fig.add_subplot(grid[4,1])
ax.set_title('CDINFB',color=mycolors[8],fontsize=title_font_size)
for i in range(10):
    ax.plot(t_ML, cAMPis_ML[0,i], alpha=0.2, color='grey',linewidth=trace_width)  
ax.plot(t_ML, np.mean(cAMPis_ML[0],axis=0), alpha=0.8, color=mycolors[8],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'Low $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylim([-0.3,1.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = fig.add_subplot(grid[5,1])
for i in range(10):
    ax.plot(t_ML, cAMPis_ML[1,i], alpha=0.2, color='grey',linewidth=trace_width)  
ax.plot(t_ML, np.mean(cAMPis_ML[1],axis=0), alpha=0.8, color=mycolors[8],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'Moderate $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.3,1.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = fig.add_subplot(grid[6,1])
for i in range(10):
    ax.plot(t_ML, cAMPis_ML[2,i], alpha=0.2, color='grey',linewidth=trace_width)  
ax.plot(t_ML, np.mean(cAMPis_ML[2],axis=0), alpha=0.8, color=mycolors[8],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'High $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.3,1.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

fig.text(0.68,0.44,'F',fontsize=letterLabelSize)
ax = fig.add_subplot(grid[4,2])
ax.set_title('IFFL',color=mycolors[0],fontsize=title_font_size)
for i in range(10):
    ax.plot(t_Kamino,ys_Kamino[0,i],alpha=0.2,color='grey',linewidth=2)
ax.plot(t_Kamino, np.mean(ys_Kamino[0],axis=0), alpha=0.8, color=mycolors[0],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'Low $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylim([-0.3,1.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = fig.add_subplot(grid[5,2])
for i in range(10):
    ax.plot(t_Kamino,ys_Kamino[1,i],alpha=0.2,color='grey',linewidth=2)
ax.plot(t_Kamino, np.mean(ys_Kamino[1],axis=0), alpha=0.8, color=mycolors[0],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'Moderate $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.3,1.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = fig.add_subplot(grid[6,2])
for i in range(10):
    ax.plot(t_Kamino,ys_Kamino[2,i],alpha=0.2,color='grey',linewidth=2)
ax.plot(t_Kamino, np.mean(ys_Kamino[2],axis=0), alpha=0.8, color=mycolors[0],linewidth=trace_width)  
ax.vlines(15,-1,10,color='k',linestyle='dashed')
ax.text(29.5,1.2,'High $cAMP_e$ Input',fontsize=12,horizontalalignment='right')
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xlim([0,30])
ax.set_xticks(np.linspace(0,30,7))
ax.set_ylim([-0.3,1.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

fig.add_artist(mpatches.FancyBboxPatch((0.033,0.04),0.93,0.41,facecolor='palegreen',zorder=0,alpha=0.5,mutation_scale=0.1,linewidth=0))

plt.subplots_adjust(top = 0.95, bottom = 0.07, right = 0.98, left = 0.055)
plt.show()


#%% Figure 4 (Population Firing Rate): pull data 

# experimental data
npzfile = np.load('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/exp_data/Gregor2010_pop_firing_rate.npz')
PopRateExp = npzfile['PopRateExp']

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Sgro_noise.pickle','rb') as f:
    firingRate_Sgro_noise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Goldbeter_noise.pickle','rb') as f:
    firingRate_Goldbeter_noise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Maeda_Noise.pickle','rb') as f:
    firingRate_Maeda_noise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Gregor_Noise.pickle','rb') as f:
    firingRate_Gregor_Noise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Kamino_Noise.pickle','rb') as f:
    firingRate_Kamino_Noise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Kamino_Noise_log.pickle','rb') as f:
    firingRate_Kamino_Noise_log = pickle.load(f)

#%% Figure 4 (Population Firing Rate): pull data 

f = plt.figure(figsize=(8,8))
gs = plt.GridSpec(3,2,wspace=0.4,hspace=0.6)

f.text(0.01,0.95,'A',fontsize=letterLabelSize)
ax = f.add_subplot(gs[0,0])
ax.set_title('Experiment',fontsize=title_font_size)
im = ax.imshow(PopRateExp,cmap='jet',aspect=1.3,vmin=0,vmax=0.18)
for i in ([2.5,3.5,4.5,5.5,6.5,7.5,8.5]):
    ax.vlines(i,-0.5,6.5,color='white')
cbar = plt.colorbar(im,ax=ax,ticks=[0,.06,.12,0.18])
cbar.set_label('cAMP Pulses/Min',fontsize=14,rotation=-90,labelpad=15)
cbar.ax.tick_params(labelsize=12)
ax.set_xlabel('Flow Rate (mL/min)', size=14)
ax.set_xticks([0,1,2,3,4,5,6,7,8,9],[0,1,2,4,6,8,10,15,16,20])
ax.set_ylabel('Cell Density (ML)', size=14)
ax.set_yticks([0,2,4,6],['1/2','1/8','1/32','1/128'],fontsize=tick_font_size)
ax.tick_params(grid_linewidth = 15, labelsize = 10)

f.add_artist(mpatches.FancyBboxPatch((0.035,0.7),0.44,0.265,facecolor='lightgrey',zorder=0,alpha=1,mutation_scale=0.1,linewidth=0))

f.text(0.52,0.95,'B',fontsize=letterLabelSize)
ax = f.add_subplot(gs[0,1])
ax.set_title('Receptor Desensitization',fontsize=title_font_size,color=mycolors[7],pad=10)
im = ax.imshow(firingRate_Goldbeter_noise,origin='lower',vmin=0,vmax=1.5,cmap='jet')
cbar = plt.colorbar(im,ax=ax)
cbar.set_label('Firing Rate (A.U.)',size=14,rotation=-90,labelpad=15)
cbar.ax.tick_params(labelsize=12)
ax.set_xlabel('Flow Rate (A.U.)',fontsize=14)
ax.set_xticks([0,12,24],[0,50,100])
ax.set_ylabel('Cell Density (A.U.)',fontsize=14)
ax.set_yticks([0,8,16,24],['$10^{-2}$','$10^{-1}$','$10^{0}$','$10^1$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)

f.text(0.01,0.62,'C',fontsize=letterLabelSize)
ax = f.add_subplot(gs[1,0])
ax.set_title('CDINFB',fontsize=title_font_size,color=mycolors[8],pad=10)
im = ax.imshow(firingRate_Maeda_noise,origin='lower',vmin=0,vmax=0.6,cmap='jet')
cbar = plt.colorbar(im,ax=ax)
cbar.set_label('Firing Rate (A.U.)',size=14,rotation=-90,labelpad=15)
cbar.ax.tick_params(labelsize=12)
ax.set_xlabel('Flow Rate (A.U.)',fontsize=14)
ax.set_xticks([0,12.5,25],[0,50,100])
ax.set_ylabel('Cell Density (A.U.)',fontsize=14)
ax.set_yticks([1,10,19,25],['$10^{0}$','$10^1$','$10^2$','$10^{2.7}$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)

f.text(0.52,0.62,'D',fontsize=letterLabelSize)
ax = f.add_subplot(gs[1,1])
ax.set_title('Phase Oscillator',fontsize=title_font_size,color=mycolors[2],pad=10)
im = ax.imshow(firingRate_Gregor_Noise,origin='lower',vmin=0,vmax=1.2,cmap='jet')
cbar = plt.colorbar(im,ax=ax,ticks=[0,.4,.8,1.2])
cbar.set_label('Firing Rate (A.U.)',size=14,rotation=-90,labelpad=15)
cbar.ax.tick_params(labelsize=12)
ax.set_xlabel('Flow Rate (A.U.)',fontsize=14)
ax.set_xticks([0,12.5,25],[0,25,50])
ax.set_ylabel('Cell Density (A.U.)',fontsize=14)
ax.set_yticks([3,9,14,19.5,25],['$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^0$','$10^1$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)

f.text(0.01,0.27,'E',fontsize=letterLabelSize)
ax = f.add_subplot(gs[2,0])
ax.set_title('IPNFB',fontsize=title_font_size,color=mycolors[5],pad=10)
im = ax.imshow(firingRate_Sgro_noise,origin='lower',vmin=0,vmax=0.6,cmap='jet')
cbar = plt.colorbar(im,ax=ax,ticks=[0,.2,.4,.6])
cbar.set_label('Firing Rate (A.U.)',size=14,rotation=-90,labelpad=15)
cbar.ax.tick_params(labelsize=12)
ax.set_xlabel('Flow Rate (A.U.)',fontsize=14)
ax.set_xticks([0,12.5,25],[0,0.5,1])
ax.set_ylabel('Cell Density (A.U.)',fontsize=14)
ax.set_yticks([0,8,17,25],['$10^{-6}$','$10^{-5}$','$10^{-4}$','$10^{-3}$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)

f.text(0.52,0.27,'F',fontsize=letterLabelSize)
ax = f.add_subplot(gs[2,1])
ax.set_title('IFFL',fontsize=title_font_size,color=mycolors[0],pad=10)
im = ax.imshow(firingRate_Kamino_Noise,origin='lower',vmin=0,vmax=0.75,cmap='jet')
cbar = plt.colorbar(im,ax=ax,ticks=[0,.25,.5,.75])
cbar.set_label('Firing Rate (A.U.)',size=14,rotation=-90,labelpad=15)
cbar.ax.tick_params(labelsize=12)
ax.set_xlabel('Flow Rate (A.U.)',fontsize=14)
ax.set_xticks([0,15,30],[0,50,100])
ax.set_ylabel('Cell Density (A.U.)',fontsize=14)
ax.set_yticks([0,15,30],['$10^{0}$','$10^{1}$','$10^{2}$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)
ax = f.add_subplot([.805,.105,.05,.05])
ax.imshow(firingRate_Kamino_Noise_log,origin='lower',vmin=0,vmax=0.75,cmap='jet')
ax.set_xticks([0,30],['$10^{0}$','$10^{2}$'],color='white')
ax.set_yticks([0,30],['$10^{0}$','$10^{2}$'])
ax.tick_params(colors='white')
ax.spines['left'].set_color('white')
ax.spines['bottom'].set_color('white')
ax.spines['right'].set_color('white')
ax.spines['top'].set_color('white')

f.add_artist(mpatches.FancyBboxPatch((0.035,0.035),0.93,0.26,facecolor='palegreen',zorder=0,alpha=0.5,mutation_scale=0.1,linewidth=0))

plt.subplots_adjust(top = 0.95, bottom = 0.07, right = 0.94, left = 0.1)
plt.show()

#%% Figure 5 (Adaptive Spiking and Bifurcation Dynamics): pull data 

# experimental data
my_dir = '//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/exp_data/'
Sgro2015Figure1excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure1')

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_adaptiveSpike_data/adaptiveSpike_Gregor.pickle','rb') as f:
   t_adapt_Gregor,cAMPi_adapt_Gregor = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_adaptiveSpike_data/adaptiveSpike_Sgro.pickle','rb') as f:
   t_adapt_Sgro,cAMPi_adapt_Sgro = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_adaptiveSpike_data/adaptiveSpike_Goldbeter.pickle','rb') as f:
   t_adapt_Goldbeter,cAMPi_adapt_Goldbeter = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_adaptiveSpike_data/adaptiveSpike_MaedaLoomis.pickle','rb') as f:
   t_adapt_Maeda,cAMPi_adapt_Maeda = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_adaptiveSpike_data/adaptiveSpike_Kamino.pickle','rb') as f:
   t_adapt_Kamino,cAMPi_adapt_Kamino = pickle.load(f)

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_bifurcationDynamics_data/bifurcationDynamics_Gregor.pickle','rb') as f:
    t_bifur_Gregor, cAMPi_bifur_Gregor = pickle.load(f) 
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_bifurcationDynamics_data/bifurcationDynamics_Sgro.pickle','rb') as f:
    t_bifur_Sgro, cAMPi_bifur_Sgro = pickle.load(f) 
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_bifurcationDynamics_data/bifurcationDynamics_Goldbeter.pickle','rb') as f:
    t_bifur_Goldbeter, cAMPi_bifur_Goldbeter = pickle.load(f) 
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_bifurcationDynamics_data/bifurcationDynamics_MaedaLoomis.pickle','rb') as f:
    t_bifur_Maeda, cAMPi_bifur_Maeda = pickle.load(f) 
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_bifurcationDynamics_data/bifurcationDynamics_Kamino.pickle','rb') as f:
    t_bifur_Kamino, cAMPi_bifur_Kamino = pickle.load(f) 

#%% Figure 5 (Adaptive Spiking and Bifurcation Dynamics): plot 

f = plt.figure(figsize=(8,8.5))

f.text(0.01,0.95,'A',fontsize=letterLabelSize)
ax = f.add_subplot([.12,.74,.35,.18])
ax.set_title('Experiment:\nAdaptive Spiking',fontsize=title_font_size)
ax.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 1 FRET Trace (1nM)"],
                              linewidth=trace_width,color='k')
ax.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 2 FRET Trace (1nM)"],
                               linewidth=trace_width,color='dimgrey')
ax.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 3 FRET Trace (1nM)"],
                               linewidth=trace_width,color='darkgrey')
ax.text(17,0.5,'1$nm$ cAMP',fontsize=tick_font_size)
ax.vlines(5,-1,1,color='k',linestyle='dashed')
ax.set_xlabel('Time (min)',fontsize=label_font_size)
ax.set_xlim([0,30])
ax.set_ylabel('FRET (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.1,.6])
ax.set_yticks(np.linspace(0,.5,2))
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.text(0.49,0.95,'B',fontsize=letterLabelSize)
ax = f.add_subplot([.55,.74,.35,.18])
ax.set_title('Experiment:\nBifurcation Dynamics',fontsize=title_font_size)
ax.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 1 FRET Trace (10uM)"],
                              linewidth=trace_width,color='k')
ax.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 2 FRET Trace (10uM)"],
                               linewidth=trace_width,color='dimgrey')
ax.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 3 FRET Trace (10uM)"],
                               linewidth=trace_width,color='darkgrey')
ax.text(15,0.5,'10$\mu m$ cAMP',fontsize=tick_font_size)
ax.vlines(5,-1,1,color='k',linestyle='dashed')
ax.set_xlabel('Time (min)',fontsize=label_font_size)
ax.set_xlim([0,30])
ax.set_ylim([-0.1,.6])
ax.set_yticks(np.linspace(0,.5,2))
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.add_artist(mpatches.FancyBboxPatch((0.033,0.69),0.935,0.275,facecolor='lightgrey',zorder=0,alpha=1,mutation_scale=0.1,linewidth=0))

f.text(0.01,0.61,'C',fontsize=letterLabelSize)
ax = f.add_subplot([.12,.4,.35,.18])
ax.set_title('Models:\nAdaptive Spiking',fontsize=title_font_size)
ax.plot(t_adapt_Goldbeter,cAMPi_adapt_Goldbeter,linewidth=2,color=mycolors[7])
ax.plot(t_adapt_Kamino,cAMPi_adapt_Kamino,linewidth=2,linestyle='dashed',color=mycolors[0])
ax.plot(t_adapt_Sgro,cAMPi_adapt_Sgro,linewidth=2,color=mycolors[5])
ax.vlines(1,-1,10,color='k',linestyle='dashed')
ax.set_xlim([0,6])
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.2,1.2])
ax.set_yticks(np.linspace(0,1,3))
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = f.add_subplot([.12,.17,.35,.18])
ax.plot(t_adapt_Maeda,cAMPi_adapt_Maeda,linewidth=2,color=mycolors[8])
ax.plot(t_adapt_Gregor,cAMPi_adapt_Gregor,linewidth=2,color=mycolors[2])
ax.vlines(1,-1,10,color='k',linestyle='dashed')
xlabel = ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
xlabel.set_position((1.15,0))
ax.set_xlim([0,6])
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.2,1.2])
ax.set_yticks(np.linspace(0,1,3))
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.text(0.49,.61,'D',fontsize=letterLabelSize)
ax = f.add_subplot([.55,.4,.35,.18])
ax.set_title('Models:\nBifurcation Dynamics',fontsize=title_font_size)
ax.plot(t_bifur_Sgro,cAMPi_bifur_Sgro,linewidth=2,color=mycolors[5])
ax.vlines(1,-1,10,color='k',linestyle='dashed')
ax.set_xlim([0,6])
ax.set_ylim([-.2,1.8])
ax.set_yticks(np.linspace(0,1.5,4))
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = f.add_subplot([.55,.17,.35,.18])
ax.plot(t_bifur_Goldbeter,cAMPi_bifur_Goldbeter,linewidth=2,color=mycolors[7])
ax.plot(t_bifur_Kamino,cAMPi_bifur_Kamino,linewidth=2,color=mycolors[0],linestyle='dashed')
ax.plot(t_bifur_Gregor,cAMPi_bifur_Gregor,linewidth=2,color=mycolors[2])
ax2 = ax.twinx()
ax2.plot(t_bifur_Maeda,cAMPi_bifur_Maeda,linewidth=2,color=mycolors[8])
ax.set_xlim([0,6])
ax.set_ylim([-.2,1.8])
ax2.set_ylim([-40,360])
ax2.set_ylabel('CDINFB $cAMP_i$ (A.U.)',fontsize=12,rotation=-90,labelpad=15)
ax.set_yticks(np.linspace(0,1.5,4))
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax2.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.add_artist(mpatches.FancyBboxPatch((0.033,0.39),0.935,0.235,facecolor='palegreen',zorder=0,alpha=0.5,mutation_scale=0.1,linewidth=0))

f.add_artist(mlines.Line2D([0.05,0.15],[0.07,0.07],color=mycolors[7],linewidth=2))
f.text(0.17,.065,'Receptor Desensitization',fontsize=12)
f.add_artist(mlines.Line2D([0.05,0.15],[0.03,0.03],color=mycolors[5],linewidth=2))
f.text(0.17,.025,'IPNFB',fontsize=12)
f.add_artist(mlines.Line2D([0.45,0.55],[0.07,0.07],color=mycolors[0],linewidth=2,linestyle='dashed'))
f.text(0.57,.065,'IFFL',fontsize=12)
f.add_artist(mlines.Line2D([0.45,0.55],[0.03,0.03],color=mycolors[8],linewidth=2))
f.text(0.57,.025,'CDINFB',fontsize=12)
f.add_artist(mlines.Line2D([0.67,0.77],[0.07,0.07],color=mycolors[2],linewidth=2))
f.text(0.79,.065,'Phase Oscillator',fontsize=12)

#%% Figure 6 (Fold Change Detection): pull data 

# experimental data
my_dir = '//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/exp_data/'
Kamino_FCD = pd.read_excel(my_dir+r'Kamino_FCD_exp_data.xlsx',sheet_name='Sheet1')

# simulation data
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_foldChange_data/FCD_sample.pickle','rb') as f:
    sample_t,sample_cAMPi = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_foldChange_data/FCD_Kamino.pickle','rb') as f:
    normPeakProm_Kamino = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_foldChange_data/FCD_Sgro.pickle','rb') as f:
    normPeakProm_Sgro = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_foldChange_data/FCD_Goldbeter.pickle','rb') as f:
    normPeakProm_Goldbeter = pickle.load(f)

FCvals = np.logspace(0,2,8)
FCvals_GB = np.logspace(0,2,12)
FCcolors = plt.cm.Greys(np.array((1,2,3,4))/4)
primingLabels = ['Priming Conc.: 0.1','Priming Conc.: 1.0','Priming Conc.: 3.0','Priming Conc.: 10.0']

#%% Figure 6 (Fold Change Detection): plot 

f = plt.figure(figsize=(8,10))
gs = plt.GridSpec(3,2,wspace=0.5,hspace=0.5)

f.text(0.01,0.96,'A',fontsize=letterLabelSize)
f.text(0.17,0.965,'Experiment: Fold-Change Detection (FCD)',fontsize=title_font_size)
ax = f.add_subplot([0.14,.86,.35,.09])
ax.plot(np.array((0,10,10,20,20,30)),np.array((0,0,1,1,2,2)),color='k',linewidth=trace_width)
ax.arrow(11,0,0,0.9,color='k',linewidth=2,length_includes_head=True,head_width=.5,head_length=.1)
ax.arrow(11,0.9,0,-0.9,color='k',linewidth=2,length_includes_head=True,head_width=.5,head_length=.1)
ax.text(12,0.1,'Priming\nConc.',fontsize=12)
ax.arrow(21,1,0,0.9,color='k',linewidth=2,length_includes_head=True,head_width=.5,head_length=.1)
ax.arrow(21,1.9,0,-0.9,color='k',linewidth=2,length_includes_head=True,head_width=.5,head_length=.1)
ax.text(22,1.1,'Fold\nChange',fontsize=12)
ax.set_xlim([0,30])
ax.set_ylabel('$cAMP_e$\nInput',fontsize=label_font_size)
ax.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax = f.add_subplot([0.14,.73,.35,.09])
ax.plot(30*sample_t/sample_t[-1],sample_cAMPi,color='k',linewidth=trace_width)
ax.arrow(22,0.08,0,.12,color='k',linewidth=2,length_includes_head=True,head_width=.5,head_length=.02)
ax.arrow(22,0.2,0,-0.12,color='k',linewidth=2,length_includes_head=True,head_width=.5,head_length=.02)
ax.text(22,0.25,'Second Peak\nProminence',fontsize=12,horizontalalignment='center')
ax.set_xlim([0,30])
ax.set_ylim([0,.4])
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_ylabel('$cAMP_i$\nResponse',fontsize=label_font_size)
ax.tick_params(axis='both', which='major', labelsize=tick_font_size)

ax = f.add_subplot(gs[0,1])
ax.errorbar(Kamino_FCD["FC_100pM"], Kamino_FCD["100pM mean"], yerr=Kamino_FCD["100pM SD"], xerr=None, color=FCcolors[0], linewidth=trace_width, label='0.1nM', ecolor=FCcolors[0], elinewidth=trace_width,capsize=5,capthick=2)
ax.errorbar(Kamino_FCD["FC_1nM"], Kamino_FCD["1nM mean"], yerr=Kamino_FCD["1nM SD"], xerr=None, color=FCcolors[1], linewidth=trace_width, label='1nM', ecolor=FCcolors[1], elinewidth=trace_width,capsize=5,capthick=2)
ax.errorbar(Kamino_FCD["FC_3nM"], Kamino_FCD["3nM mean"], yerr=Kamino_FCD["3nM SD"], xerr=None, color=FCcolors[2], linewidth=trace_width, label='3nM', ecolor=FCcolors[2], elinewidth=trace_width,capsize=5,capthick=2)
ax.errorbar(Kamino_FCD["FC_10nM"], Kamino_FCD["10nM mean"], yerr=Kamino_FCD["10nM SD"], xerr=None, color=FCcolors[3], linewidth=trace_width, label='10nM', ecolor=FCcolors[3], elinewidth=trace_width,capsize=5,capthick=2)
ax.legend(frameon=0,ncol=2,fontsize=12)
ax.set_ylim([1.5,2.6])
ax.set_ylabel('Second Peak\nProminence (A.U.)',fontsize=label_font_size)
ax.set_xlabel(r'$cAMP_{e}$'+' fold change',fontsize=label_font_size)
ax.set_xscale('log')
ax.tick_params(axis='both', which='major', labelsize=tick_font_size)

f.add_artist(mpatches.FancyBboxPatch((0.033,0.69),0.935,0.277,facecolor='lightgrey',zorder=0,alpha=1,mutation_scale=0.1))

f.text(0.01,0.61,'B',fontsize=letterLabelSize)
ax = f.add_subplot(gs[1,0])
ax.set_title('IFFL',color=mycolors[0],fontsize=title_font_size)
for i in range(len(normPeakProm_Kamino)):
    ax.plot(FCvals,normPeakProm_Kamino[i],color=plt.cm.Greys(i/len(normPeakProm_Kamino)),linewidth=trace_width)
ax.set_xscale('log')
ax.set_xlabel('$cAMP_e$ Fold Change',fontsize=label_font_size)
ax.set_ylabel('Second Peak\nProminence (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.1,1.1])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.text(0.51,0.61,'C',fontsize=letterLabelSize)
ax = f.add_subplot(gs[1,1])
ax.set_title('IPNFB',color=mycolors[5],fontsize=title_font_size)
for i in range(len(normPeakProm_Sgro)):
    ax.plot(FCvals,np.mean(normPeakProm_Sgro[i],axis=1),color=plt.cm.Greys(i/len(normPeakProm_Kamino)),linewidth=trace_width)
    ax.errorbar(FCvals,np.mean(normPeakProm_Sgro[i],axis=1),np.std(normPeakProm_Sgro[i],axis=1)/np.sqrt(np.shape(normPeakProm_Sgro)[2]),capsize=5,capthick=2,color=plt.cm.Greys(i/len(normPeakProm_Kamino)),linewidth=trace_width)
ax.set_xscale('log')
ax.set_xlabel('$cAMP_e$ Fold Change',fontsize=label_font_size)
ax.set_ylabel('Second Peak\nProminence (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.1,1.1])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.add_artist(mpatches.FancyBboxPatch((0.033,0.36),0.935,0.265,facecolor='palegreen',zorder=0,alpha=0.5,mutation_scale=0.1))


f.text(0.01,0.29,'D',fontsize=letterLabelSize)
ax = f.add_subplot(gs[2,0])
ax.set_title('Receptor Desensitization',color=mycolors[7],fontsize=title_font_size)
for i in range(len(normPeakProm_Goldbeter)):
    ax.plot(FCvals_GB,normPeakProm_Goldbeter[i],color=plt.cm.Greys((i+1)/(len(normPeakProm_Goldbeter)+1)),linewidth=trace_width,label=primingLabels[i])
ax.set_xscale('log')
ax.set_xlabel('$cAMP_e$ Fold Change',fontsize=label_font_size)
ax.set_ylabel('Second Peak\nProminence (A.U.)',fontsize=label_font_size)
ax.set_ylim([-1,51.1])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax.legend(frameon=0,fontsize=tick_font_size,loc='upper left',bbox_to_anchor=[1.3,.4,.5,.5])

plt.subplots_adjust(top = 0.95, bottom = 0.07, right = 0.97, left = 0.15)
plt.show()

#%% Figure S1 (Population Firing Rate, noise and no noise): pull data 

# experimental data
npzfile = np.load('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/exp_data/Gregor2010_pop_firing_rate.npz')
PopRateExp = npzfile['PopRateExp']


with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Sgro_noise.pickle','rb') as f:
    firingRate_Sgro_noise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Goldbeter_noise.pickle','rb') as f:
    firingRate_Goldbeter_noise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Maeda_Noise.pickle','rb') as f:
    firingRate_Maeda_noise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Gregor_Noise.pickle','rb') as f:
    firingRate_Gregor_Noise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Kamino_Noise.pickle','rb') as f:
    firingRate_Kamino_Noise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Kamino_Noise_log.pickle','rb') as f:
    firingRate_Kamino_Noise_log = pickle.load(f)

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Sgro_noNoise.pickle','rb') as f:
    firingRate_Sgro_noNoise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Goldbeter_noNoise.pickle','rb') as f:
    firingRate_Goldbeter_noNoise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Maeda_noNoise.pickle','rb') as f:
    firingRate_Maeda_noNoise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Gregor_noNoise.pickle','rb') as f:
    firingRate_Gregor_noNoise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Kamino_noNoise.pickle','rb') as f:
    firingRate_Kamino_noNoise = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Kamino_noNoise_log.pickle','rb') as f:
    firingRate_Kamino_noNoise_log = pickle.load(f)

#%% Figure S1

f = plt.figure(figsize=(16,9))
gs = plt.GridSpec(3,5,wspace=0.5,hspace=0.3)

f.text(0.01,0.95,'A',fontsize=letterLabelSize)
ax = f.add_subplot(gs[0,0:2])
ax.set_title('Experiment',fontsize=label_font_size)
im = ax.imshow(PopRateExp,cmap='jet',aspect=0.6,vmin=0,vmax=0.18)
for i in ([2.5,3.5,4.5,5.5,6.5,7.5,8.5]):
    ax.vlines(i,-0.5,6.5,color='white')
cbar = plt.colorbar(im,ax=ax,ticks=[0,.06,.12,0.18],shrink=0.75,aspect=10)
cbar.set_label('cAMP Pulses/Min',fontsize=14,rotation=-90,labelpad=15)
cbar.ax.tick_params(labelsize=12)
ax.set_xlabel('Flow Rate (mL/min)', size=14)
ax.set_xticks([0,1,2,3,4,5,6,7,8,9],[0,1,2,4,6,8,10,15,16,20])
ax.set_ylabel('Cell Density (ML)', size=14)
ax.set_yticks([0,2,4,6],['1/2','1/8','1/32','1/128'],fontsize=tick_font_size)
ax.tick_params(grid_linewidth = 15, labelsize = 12)

f.add_artist(mpatches.FancyBboxPatch((0.02,0.71),0.37,0.27,facecolor='lightgrey',zorder=0,alpha=1,mutation_scale=0.05,linewidth=0))

f.text(0.01,.65,'B',fontsize=letterLabelSize)
ax = f.add_subplot(gs[1,0])
ax.set_title('Receptor Desensitization\nwith Noise',fontsize=label_font_size,color=mycolors[7],pad=10)
im = ax.imshow(firingRate_Goldbeter_noise,origin='lower',vmin=0,vmax=1.5,cmap='jet')
cbar = plt.colorbar(im,ax=ax,shrink=0.72,aspect=10)
cbar.ax.tick_params(labelsize=12)
ax.set_xticks([0,12,24],[0,50,100])
ax.set_yticks([0,8,16,24],['$10^{-2}$','$10^{-1}$','$10^{0}$','$10^1$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)

ax = f.add_subplot(gs[2,0])
ax.set_title('Receptor Desensitization\nwithout Noise',fontsize=label_font_size,color=mycolors[7],pad=10)
im = ax.imshow(firingRate_Goldbeter_noNoise,origin='lower',vmin=0,vmax=1.5,cmap='jet')
cbar = plt.colorbar(im,ax=ax,shrink=0.72,aspect=10)
cbar.ax.tick_params(labelsize=12)
ax.set_xticks([0,12,24],[0,50,100])
ax.set_yticks([0,8,16,24],['$10^{-2}$','$10^{-1}$','$10^{0}$','$10^1$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)

f.text(0.23,.65,'C',fontsize=letterLabelSize)
ax = f.add_subplot(gs[1,1])
ax.set_title('CDINFB\nwith Noise',fontsize=label_font_size,color=mycolors[8],pad=10)
im = ax.imshow(firingRate_Maeda_noise,origin='lower',vmin=0,vmax=0.6,cmap='jet')
cbar = plt.colorbar(im,ax=ax,shrink=0.72,aspect=10)
cbar.ax.tick_params(labelsize=12)
ax.set_xticks([0,12.5,25],[0,50,100])
ax.set_yticks([1,10,19,25],['$10^{0}$','$10^1$','$10^2$','$10^{2.7}$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)

ax = f.add_subplot(gs[2,1])
ax.set_title('CDINFB\nwithout Noise',fontsize=label_font_size,color=mycolors[8],pad=10)
im = ax.imshow(firingRate_Maeda_noNoise,origin='lower',vmin=0,vmax=0.6,cmap='jet')
cbar = plt.colorbar(im,ax=ax,shrink=0.72,aspect=10)
cbar.ax.tick_params(labelsize=12)
ax.set_xticks([0,12.5,25],[0,50,100])
ax.set_yticks([1,10,19,25],['$10^{0}$','$10^1$','$10^2$','$10^{2.7}$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)

f.text(0.42,.65,'D',fontsize=letterLabelSize)
ax = f.add_subplot(gs[1,2])
ax.set_title('Phase Oscillator\nwith Noise',fontsize=label_font_size,color=mycolors[2],pad=10)
im = ax.imshow(firingRate_Gregor_Noise,origin='lower',vmin=0,vmax=1.2,cmap='jet')
cbar = plt.colorbar(im,ax=ax,ticks=[0,.4,.8,1.2],shrink=0.72,aspect=10)
cbar.ax.tick_params(labelsize=12)
ax.set_xticks([0,12.5,25],[0,25,50])
ax.set_yticks([3,9,14,19.5,25],['$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^0$','$10^1$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)

ax = f.add_subplot(gs[2,2])
ax.set_title('Phase Oscillator\nwithout Noise',fontsize=label_font_size,color=mycolors[2],pad=10)
im = ax.imshow(firingRate_Gregor_noNoise,origin='lower',vmin=0,vmax=1.2,cmap='jet')
cbar = plt.colorbar(im,ax=ax,ticks=[0,.4,.8,1.2],shrink=0.72,aspect=10)
cbar.ax.tick_params(labelsize=12)
ax.set_xticks([0,12.5,25],[0,25,50])
ax.set_yticks([3,9,14,19.5,25],['$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^0$','$10^1$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)

f.text(0.61,.65,'E',fontsize=letterLabelSize)
ax = f.add_subplot(gs[1,3])
ax.set_title('IPNFB\nwith Noise',fontsize=label_font_size,color=mycolors[5],pad=10)
im = ax.imshow(firingRate_Sgro_noise,origin='lower',vmin=0,vmax=0.6,cmap='jet')
cbar = plt.colorbar(im,ax=ax,ticks=[0,.2,.4,.6],shrink=0.72,aspect=10)
cbar.ax.tick_params(labelsize=12)
ax.set_xticks([0,12.5,25],[0,0.5,1])
ax.set_yticks([0,8,17,25],['$10^{-6}$','$10^{-5}$','$10^{-4}$','$10^{-3}$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)

ax = f.add_subplot(gs[2,3])
ax.set_title('IPNFB\nwithout Noise',fontsize=label_font_size,color=mycolors[5],pad=10)
im = ax.imshow(firingRate_Sgro_noNoise,origin='lower',vmin=0,vmax=0.6,cmap='jet')
cbar = plt.colorbar(im,ax=ax,ticks=[0,.2,.4,.6],shrink=0.72,aspect=10)
cbar.ax.tick_params(labelsize=12)
ax.set_xticks([0,12.5,25],[0,0.5,1])
ax.set_yticks([0,8,17,25],['$10^{-6}$','$10^{-5}$','$10^{-4}$','$10^{-3}$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)

f.text(0.81,0.65,'F',fontsize=letterLabelSize)
ax = f.add_subplot(gs[1,4])
ax.set_title('IFFL\nwith Noise',fontsize=title_font_size,color=mycolors[0],pad=10)
im = ax.imshow(firingRate_Kamino_Noise,origin='lower',vmin=0,vmax=0.75,cmap='jet')
cbar = plt.colorbar(im,ax=ax,ticks=[0,.25,.5,.75],shrink=0.72,aspect=10)
cbar.ax.tick_params(labelsize=12)
ax.set_xticks([0,15,30],[0,50,100])
ax.set_yticks([0,15,30],['$10^{0}$','$10^{1}$','$10^{2}$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)
ax = f.add_subplot([.898,.465,.05,.05])
ax.imshow(firingRate_Kamino_Noise_log,origin='lower',vmin=0,vmax=0.75,cmap='jet')
ax.set_xticks([0,30],['$10^{0}$','$10^{2}$'],color='white')
ax.set_yticks([0,30],['$10^{0}$','$10^{2}$'])
ax.tick_params(colors='white')
ax.spines['left'].set_color('white')
ax.spines['bottom'].set_color('white')
ax.spines['right'].set_color('white')
ax.spines['top'].set_color('white')

ax = f.add_subplot(gs[2,4])
ax.set_title('IFFL\nwithout Noise',fontsize=title_font_size,color=mycolors[0],pad=10)
im = ax.imshow(firingRate_Kamino_noNoise,origin='lower',vmin=0,vmax=0.6,cmap='jet')
cbar = plt.colorbar(im,ax=ax,ticks=[0,.25,.5,.75],shrink=0.72,aspect=10)
cbar.ax.tick_params(labelsize=12)
ax.set_xticks([0,15,30],[0,50,100])
ax.set_yticks([0,15,30],['$10^{0}$','$10^{1}$','$10^{2}$'])
ax.tick_params(grid_linewidth = 15, labelsize = 12)
ax = f.add_subplot([.898,.138,.05,.05])
ax.imshow(firingRate_Kamino_noNoise_log,origin='lower',vmin=0,vmax=0.75,cmap='jet')
ax.set_xticks([0,30],['$10^{0}$','$10^{2}$'],color='white')
ax.set_yticks([0,30],['$10^{0}$','$10^{2}$'])
ax.tick_params(colors='white')
ax.spines['left'].set_color('white')
ax.spines['bottom'].set_color('white')
ax.spines['right'].set_color('white')
ax.spines['top'].set_color('white')

f.text(0.02,0.2,'Population Density (A.U.)',fontsize=title_font_size,rotation=90)
f.text(0.46,0.04,'Flow Rate (A.U.)',fontsize=title_font_size)

plt.subplots_adjust(top = 0.98, bottom = 0.07, right = 0.97, left = 0.08)
plt.show()

#%% Figure S2 (Sgro Model, FCD): pull data 

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_foldChange_data/FCD_Sgro.pickle','rb') as f:
    normPeakProm_Sgro = pickle.load(f)
    
#%% Figure S2 (Sgro Model, FCD): plot 

f = plt.figure(figsize=(16,4))
gs = plt.GridSpec(1,4,wspace=0.3,hspace=0.5)

f.text(0.49,0.9,'IPNFB',color=mycolors[5],fontsize=letterLabelSize)
f.text(0.46,0.02,'$cAMP_e$ Fold Change',fontsize=label_font_size)

for i in range(len(normPeakProm_Sgro)):
    ax = f.add_subplot(gs[0,i])
    ax.set_title('Priming Conc.: 0.01',fontsize=title_font_size)
    parts = ax.violinplot(np.transpose(normPeakProm_Sgro[i]),showmeans=True,showextrema=False)
    for pc in parts['bodies']:
        pc.set_facecolor(FCcolors[i])
    ax.set_xticks(np.linspace(1,8,8),np.round(FCvals).astype('int'))
    ax.set_ylim([0,1.1])
    if i==0: ax.set_ylabel('Second Peak\nProminence (A.U.)',fontsize=label_font_size)
    ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

plt.subplots_adjust(top = 0.8, bottom = 0.16, right = 0.99, left = 0.08)
plt.show()

#%% Figure S3 (Step vs Ramp): pull data

my_dir = '//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/exp_data/'
Sgro2015Figure3excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure3')

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_stepVsRamp_data/stepVsRamp_Goldbeter.pickle','rb') as f:
    t_Goldbeter,cAMP_Goldbeter = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_stepVsRamp_data/stepVsRamp_Sgro.pickle','rb') as f:
    t_Sgro,cAMP_Sgro = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_stepVsRamp_data/stepVsRamp_MaedaLoomis.pickle','rb') as f:
    t_Maeda,cAMP_Maeda = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_stepVsRamp_data/stepVsRamp_Kamino.pickle','rb') as f:
    t_Kamino,cAMP_Kamino = pickle.load(f)


#%% Figure S3 (Step vs Ramp): plot 

f = plt.figure(figsize=(16,3.6))
gs = plt.GridSpec(3,5,wspace=0.4,hspace=0.1)

f.text(0.001,0.86,'A',fontsize=letterLabelSize)
ax = f.add_subplot(gs[0, 0])
ax.plot(Sgro2015Figure3excel["Ramp Input (min Time)"],Sgro2015Figure3excel["Ramp Input (nM cAMP)"],
                              color='k', linewidth = trace_width)
ax.set_title('Experiment',color='k', fontsize=title_font_size)
ax.set_ylabel('$cAMP_e$\n(nM)',fontsize=label_font_size,labelpad=0)
ax.set_xlim([0,80])
ax.set_xticks([0,20,40,60,80],[])
ax.set_ylim([-0.1,1.1])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = f.add_subplot(gs[1:3,0])
ax.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 1 FRET Trace"],
                                color='k', linewidth = trace_width)
ax.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 2 FRET Trace"],
                                color='grey', linewidth = trace_width)
ax.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 3 FRET Trace"],
                                color='lightgrey', linewidth = trace_width)
ax.axvspan(10,30,facecolor='grey',alpha=0.2)
ax.axvspan(50,70,facecolor='grey',alpha=0.2)
ax.set_xlabel('Time (min)',fontsize=label_font_size)
ax.set_xlim([0,80])
ax.set_xticks([0,20,40,60,80])
ax.set_ylabel('FRET (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.1,.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.add_artist(mpatches.FancyBboxPatch((0.032,0.032),0.14,0.93,facecolor='lightgrey',zorder=0,alpha=1,mutation_scale=0.1,linewidth=0))

f.text(0.2,0.86,'B',fontsize=letterLabelSize)
ax = f.add_subplot(gs[0:3,1])
ax.set_title('Receptor\nDesensitization',fontsize=title_font_size,color=mycolors[7])
ax.plot(t_Goldbeter,cAMP_Goldbeter,linewidth=trace_width,color=mycolors[7])
ax.axvspan(2,6,facecolor='grey',alpha=0.2)
ax.axvspan(10,14,facecolor='grey',alpha=0.2)
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xticks([0,4,8,12,16])
ax.set_xlim([0,16])
ax.set_ylabel('$cAMP_i$ (A.U.)',fontsize=label_font_size)
ax.set_ylim([-0.2,1.1])
ax.set_yticks(np.linspace(0,1,6))
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.text(0.41,0.86,'C',fontsize=letterLabelSize)
ax = f.add_subplot(gs[0:3,2])
ax.set_title('IPNFB',fontsize=title_font_size,color=mycolors[5])
ax.plot(t_Sgro,cAMP_Sgro,linewidth=trace_width,color=mycolors[5])
ax.axvspan(2,6,facecolor='grey',alpha=0.2)
ax.axvspan(10,14,facecolor='grey',alpha=0.2)
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xticks([0,4,8,12,16])
ax.set_xlim([0,16])
ax.set_ylim([-0.2,1.1])
ax.set_yticks(np.linspace(0,1,6))
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.add_artist(mpatches.FancyBboxPatch((0.435,0.032),0.14,0.93,facecolor='lightgreen',zorder=0,alpha=0.25,mutation_scale=0.1,linewidth=0))

f.text(0.61,0.86,'D',fontsize=letterLabelSize)
ax = f.add_subplot(gs[0:3,3])
ax.set_title('CDINFB',fontsize=title_font_size,color=mycolors[8])
ax.plot(t_Maeda,cAMP_Maeda,linewidth=trace_width,color=mycolors[8])
ax.axvspan(2,6,facecolor='grey',alpha=0.2)
ax.axvspan(10,14,facecolor='grey',alpha=0.2)
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xticks([0,4,8,12,16])
ax.set_xlim([0,16])
ax.set_ylim([-0.2,1.1])
ax.set_yticks(np.linspace(0,1,6))
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.text(0.81,0.86,'E',fontsize=letterLabelSize)
ax = f.add_subplot(gs[0:3,4])
ax.set_title('IFFL',fontsize=title_font_size,color=mycolors[0])
ax.plot(t_Maeda,cAMP_Maeda,linewidth=trace_width,color=mycolors[0])
ax.axvspan(2,6,facecolor='grey',alpha=0.2)
ax.axvspan(10,14,facecolor='grey',alpha=0.2)
ax.set_xlabel('Time (A.U.)',fontsize=label_font_size)
ax.set_xticks([0,4,8,12,16])
ax.set_xlim([0,16])
ax.set_ylim([-0.2,1.1])
ax.set_yticks(np.linspace(0,1,6))
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

plt.subplots_adjust(top = 0.84, bottom = 0.2, right = 0.99, left = 0.05)
plt.show()

#%% Figure S4 (Entrainment): pull data 

my_dir = '//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/exp_data/'
Sgro2015Figure3excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure3')
Sgro2015Figure4excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure4traces')
entrainmentRs = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure4heatmap')
entrainmentRs = entrainmentRs.to_numpy()
entrainmentRs[np.where(entrainmentRs==0)] = 'NaN'

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_Entrainment_data/entrainment_Goldbeter.pickle','rb') as f:
    Rs_Goldbeter = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_Entrainment_data/entrainment_Sgro.pickle','rb') as f:
    Rs_Sgro = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_Entrainment_data/entrainment_Maeda.pickle','rb') as f:
    Rs_Maeda = pickle.load(f)
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_Entrainment_data/entrainment_Kamino.pickle','rb') as f:
    Rs_Kamino = pickle.load(f)

#%% Figure S4 (Entrainment): plot 

f = plt.figure(figsize=(8,10))
gs = plt.GridSpec(3,2,wspace=0.4,hspace=0.5)

f.text(0.01,0.95,'A',fontsize=letterLabelSize)
f.text(0.4,0.97,'Experiment',fontsize=title_font_size)
ax = f.add_subplot([0.11,.86,.4,.07])
ax.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 1 FRET Trace (1 min pulse)"],color='k',linewidth=trace_width)
ax.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 2 FRET Trace (1 min pulse)"],color='grey',linewidth=trace_width)
ax.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 3 FRET Trace (1 min pulse)"],color='lightgrey',linewidth=trace_width)
for i in range(3):
    ax.axvspan(5+i*6, 5+i*6+1, alpha=0.2, color='grey')
ax.set_xlim([0,23.5])
ylabel = ax.set_ylabel('FRET (A.U.)',fontsize=label_font_size)
ylabel.set_position((.1,-.2))
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = f.add_subplot([0.11,.76,.4,.07])
ax.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 1 FRET Trace (5 min pulse)"],color='k',linewidth=trace_width)
ax.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 2 FRET Trace (5 min pulse)"],color='grey',linewidth=trace_width)
ax.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 3 FRET Trace (5 min pulse)"],color='lightgrey',linewidth=trace_width)
for i in range(3): ax.axvspan(6+i*6, 6+i*6+5, alpha=0.2, color='grey')
ax.set_xlim([0,23.4])
ax.set_xlabel('Time (min)',fontsize=label_font_size)
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

ax = f.add_subplot(gs[0,1])
im = ax.imshow(entrainmentRs,vmin=0,vmax=1,cmap='jet')
cbar = plt.colorbar(im,ax=ax,ticks=[0,.5,1])
cbar.set_label('Entrainment Quality',size=label_font_size,rotation=-90,labelpad=20)
cbar.ax.tick_params(labelsize=tick_font_size)
ax.set_xlabel('Period (min)',fontsize=label_font_size)
ax.set_xticks([0,1,2,3],[3,4,5,6])
ax.set_ylabel('Peak Width (min)',fontsize=label_font_size)
ax.set_yticks([0,1,2,3],[4,3,2,1])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.add_artist(mpatches.FancyBboxPatch((0.032,0.69),0.935,0.275,facecolor='lightgrey',zorder=0,alpha=1,mutation_scale=0.1,linewidth=0))

f.text(0.01,.61,'B',fontsize=letterLabelSize)
ax = f.add_subplot(gs[1,0])
ax.set_title('Receptor Desensitization',fontsize=title_font_size,color=mycolors[7],pad=10)
im = ax.imshow(Rs_Goldbeter,vmin=0,vmax=1,cmap='jet',origin='lower',aspect=0.8)
cbar = plt.colorbar(im,ax=ax,ticks=[0,.5,1])
cbar.ax.tick_params(labelsize=tick_font_size)
ax.set_xticks([0,7],[0.8,1.8])
ax.set_xlim([-0.5,7.5])
ax.set_yticks([0,9],[0.3,1.4])
ax.set_ylim([-0.5,9.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.text(0.07,.3,'C',fontsize=letterLabelSize)
ax = f.add_subplot(gs[2,0])
ax.set_title('IPNFB',fontsize=title_font_size,color=mycolors[5])
im = ax.imshow(Rs_Sgro,vmin=0,vmax=1,cmap='jet',origin='lower',aspect=0.8)
cbar = plt.colorbar(im,ax=ax,ticks=[0,.5,1])
cbar.ax.tick_params(labelsize=tick_font_size)
xlabel = ax.set_xlabel('Period (A.U.)',fontsize=label_font_size)
xlabel.set_position((1.5,0))
ylabel = ax.set_ylabel('Peak Width (A.U.)',fontsize=label_font_size)
ylabel.set_position((0,1.3))
ax.set_xticks([0,7],[0.8,1.8])
ax.set_xlim([-0.5,7.5])
ax.set_yticks([0,9],[0.3,1.4])
ax.set_ylim([-0.5,9.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.text(0.51,.61,'D',fontsize=letterLabelSize)
ax = f.add_subplot(gs[1,1])
ax.set_title('CDINFB',fontsize=title_font_size,color=mycolors[8])
im = ax.imshow(Rs_Maeda,vmin=0,vmax=1,cmap='jet',origin='lower',aspect=0.8)
cbar = plt.colorbar(im,ax=ax,ticks=[0,.5,1])
cbar.ax.tick_params(labelsize=tick_font_size)
ax.set_xticks([0,7],[0.8,1.8])
ax.set_xlim([-0.5,7.5])
ax.set_yticks([0,9],[0.3,1.4])
ax.set_ylim([-0.5,9.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.text(0.51,.3,'E',fontsize=letterLabelSize)
ax = f.add_subplot(gs[2,1])
ax.set_title('IFFL',fontsize=title_font_size,color=mycolors[0])
im = ax.imshow(Rs_Kamino,vmin=0,vmax=1,cmap='jet',origin='lower',aspect=0.8)
cbar = plt.colorbar(im,ax=ax,ticks=[0,.5,1])
cbar.ax.tick_params(labelsize=tick_font_size)
ax.set_xticks([0,7],[0.8,1.8])
ax.set_xlim([-0.5,7.5])
ax.set_yticks([0,9],[0.3,1.4])
ax.set_ylim([-0.5,9.5])
ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)

f.add_artist(mpatches.FancyBboxPatch((0.032,0.032),0.935,0.595,facecolor='lightgreen',zorder=0,alpha=0.25,mutation_scale=0.1,linewidth=0))

plt.subplots_adjust(top = 0.95, bottom = 0.07, right = 0.94, left = 0.13)
plt.show()