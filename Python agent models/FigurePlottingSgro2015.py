# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 15:31:28 2019

@author: Xinwen
"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd

#Sgro2015Figure1 = scipy.io.loadmat('Sgro2015dataandfigures/DataFigure1.mat')
#Sgro2015Figure3 = scipy.io.loadmat('Sgro2015dataandfigures/Figure3Data.mat')
Sgro2015Figure4 = scipy.io.loadmat('Sgro2015dataandfigures/Figure4Data.mat')
Sgro2015Figure5 = scipy.io.loadmat('Sgro2015dataandfigures/Fig5Data.mat')
#Sgro2015Figure6 = scipy.io.loadmat('Sgro2015dataandfigures/Figure6Data.mat')
Sgro2015Figure1excel = pd.read_excel(r'Sgro2015dataandfigures/Sgro2015DataFormattedforPython.xlsx',sheetname='Figure1')
Sgro2015Figure3excel = pd.read_excel(r'Sgro2015dataandfigures/Sgro2015DataFormattedforPython.xlsx',sheetname='Figure3')
Sgro2015Figure4excel = pd.read_excel(r'Sgro2015dataandfigures/Sgro2015DataFormattedforPython.xlsx',sheetname='Figure4')
Sgro2015Figure6excel = pd.read_excel(r'Sgro2015dataandfigures/Sgro2015DataFormattedforPython.xlsx',sheetname='Figure6')


#%% Figure 1
title_font_size = 20
label_font_size = 20
tick_font_size = 16
legend_font_size = 12
trace_width = 5


fig1 = plt.figure(figsize=(12,13))
grid = plt.GridSpec(1, 2, wspace=0.2, hspace=0.2)
ax1= fig1.add_subplot(grid[0, 0])

ax1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 1 FRET Trace (1nM)"])
ax1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 2 FRET Trace (1nM)"])
ax1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 3 FRET Trace (1nM)"])

ax1.axvline(x=5, ls='--') #dashed line at 5 (cAMP onset)

# ax1.set_ylim([-0.3,1])
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('1nM cAMP', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})


ax2= fig1.add_subplot(grid[0, 1])

ax2.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 1 FRET Trace (10uM)"])
ax2.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 2 FRET Trace (10uM)"])
ax2.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 3 FRET Trace (10uM)"])

ax2.axvline(x=5, ls='--')

ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('10 uM cAMP', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
plt.show()


#%% Figure 3
fig2 = plt.figure(figsize=(12,13))
grid = plt.GridSpec(2, 1, wspace=0.2, hspace=0.2)
ax1= fig2.add_subplot(grid[0, 0])

#need to generate the times?
#Fig3Atimes = [x*80/Sgro2015Figure3["ramp1nMinput"].size for x in range(0,Sgro2015Figure3["ramp1nMinput"].size)]
#Fig3Atimes_output = [x*80/Sgro2015Figure3["ramp1nMoutput1"].size for x in range(0,Sgro2015Figure3["ramp1nMoutput1"].size)]
#ax1.plot(Fig3Atimes,Sgro2015Figure3["ramp1nMinput"][0,:])
ax1.plot(Sgro2015Figure3excel["Ramp Input (min Time)"],Sgro2015Figure3excel["Ramp Input (nM cAMP)"])

ax2= fig2.add_subplot(grid[1, 0])
ax2.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 1 FRET Trace"])
ax2.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 2 FRET Trace"])
ax2.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 3 FRET Trace"])
#ax2.plot(Fig3Atimes_output,Sgro2015Figure3["ramp1nMoutput1"][0,:])
#ax2.plot(Fig3Atimes_output,Sgro2015Figure3["ramp1nMoutput2"][0,:])
#ax2.plot(Fig3Atimes_output,Sgro2015Figure3["ramp1nMoutput3"][0,:])
#%% Figure 4
fig3 = plt.figure(figsize=(12,13))
grid = plt.GridSpec(3, 1, wspace=0.2, hspace=0.2)
ax1= fig3.add_subplot(grid[0, 0])

#need to generate the times?
#Fig4Atimes_good = [x*60/271 for x in range(0,271)]
#Fig4Atimes_bad = [x*60/271 for x in range(0,271)]

#ax1.plot(Fig4Atimes_good[0:160],Sgro2015Figure4["goodentrainment"][0,0:160])
#ax1.plot(Fig4Atimes_good[0:160],Sgro2015Figure4["goodentrainment"][1,0:160])
#ax1.plot(Fig4Atimes_good[0:160],Sgro2015Figure4["goodentrainment"][2,0:160])
ax1.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 1 FRET Trace (1 min pulse)"])
ax1.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 2 FRET Trace (1 min pulse)"])
ax1.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 3 FRET Trace (1 min pulse)"])
ax1.set_xlim([0,30])

ax2= fig3.add_subplot(grid[1, 0])
#ax2.plot(Fig4Atimes_bad[0:160],Sgro2015Figure4["badentrainment"][0,0:160])
#ax2.plot(Fig4Atimes_bad[0:160],Sgro2015Figure4["badentrainment"][1,0:160])
#ax2.plot(Fig4Atimes_bad[0:160],Sgro2015Figure4["badentrainment"][2,0:160])
ax2.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 1 FRET Trace (5 min pulse)"])
ax2.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 2 FRET Trace (5 min pulse)"])
ax2.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 3 FRET Trace (5 min pulse)"])
ax2.set_xlim([0,30])

#heatmap of Sgro2015Figure4["entrainmentRs"]
ax3= fig3.add_subplot(grid[2, 0])
ax3=plt.imshow(Sgro2015Figure4["entrainmentRs"][:,0:4], interpolation='nearest')

#%% Figure 5
fig4 = plt.figure(figsize=(12,13)) #need data from Gregor2010
plt.imshow(Sgro2015Figure5["cAMP_PulsesPerMin"], interpolation='nearest')

#%% Figure 6
fig5 = plt.figure(figsize=(12,13))
grid = plt.GridSpec(5, 1, wspace=0.2, hspace=0.2)


#void384 data = MATLAB structure, need to figure out what is what

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["No External cAMP Mean Trace"])
ax2= fig5.add_subplot(grid[1, 0])
ax2.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["Low External cAMP Mean Trace"])
ax3= fig5.add_subplot(grid[2, 0])
ax3.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["Intermediate External cAMP Mean Trace"])
ax4= fig5.add_subplot(grid[3, 0])
ax4.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["High External cAMP Mean Trace"])



#%% Figure 6 just sustained oscillation
fig6 = plt.figure(figsize=(12,13))
grid = plt.GridSpec(5, 1, wspace=0.2, hspace=0.2)
ax1= fig6.add_subplot(grid[0, 0])
ax1.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["No External cAMP Mean Trace"])
ax1.set_xlim([0,60])