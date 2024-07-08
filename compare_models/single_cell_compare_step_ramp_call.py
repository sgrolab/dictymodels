# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:23:10 2019

@author: Chuqiao
"""

import numpy as np
import random
import math
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd


from NormParam import *

##%% Get ramp input from experimental data
#Exp_time = Sgro2015Figure3excel["Ramp Input (min Time)"]# np.linspace(0,20,num=1000) # t step 0.02
#Exp_time = Exp_time[~np.isnan(Exp_time)]
#RampInput_Exp = Sgro2015Figure3excel["Ramp Input (nM cAMP)"]
#RampInput_Exp = RampInput_Exp[~np.isnan(RampInput_Exp)]

# experimental data
my_dir = '//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/exp_data/'
Sgro2015Figure3excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure3')


# set matplotlib default font
import matplotlib
font = {'family' : 'Roboto'}
matplotlib.rc('font', **font)

#%% Single cell step cs. ramp inputs
# Sgro 2015
from Sgro2015_agent_and_pop_FUN import Sgro2015_agent

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
Nt=27; # normalization factor of t
Na=3.5;  # normalization factor of A
A0=-1.5; R0=-0.5
Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)

## stimulation settings
#cAMP=1
T_tot = 16; T_tot=T_tot*Nt_Sgro
dt=0.001 ;  t=list(np.arange(0,T_tot,dt))
#signal_trace_Sgro=np.zeros(len(t))
#
## define step stimulus
#TStepOn = 2.5; tStepOn = math.floor(TStepOn/dt*Nt_Sgro); 
#TStepOff = 7.5; tStepOff = math.floor(TStepOff/dt*Nt_Sgro); 
#signal_trace_Sgro[tStepOn:tStepOff] = cAMP
#
## define exponential ramp stimulus
##TRampOn=10; tRampOn=math.floor(TRampOn/dt*Nt_Sgro)
##TRampOff=24; tRampOff=math.floor(TRampOff/dt*Nt_Sgro)
#TRampOn=12.5; tRampOn=math.floor(TRampOn/dt*Nt_Sgro)
#TRampOff=17.5; tRampOff=math.floor(TRampOff/dt*Nt_Sgro)######################## 10/23 here!!!!
#tRamp=np.linspace(tRampOn, tRampOff, num=tRampOff-tRampOn+1)
#
#a=math.log(cAMP+1)
#temp_t=np.linspace(0,1,len(tRamp))
#signal_trace_Sgro[tRampOn:tRampOff+1] = np.exp(a*temp_t)-1
cAMP = 1
Exp_time = Sgro2015Figure3excel["Ramp Input (min Time)"]# np.linspace(0,20,num=1000) # t step 0.02
Exp_time = Exp_time[~np.isnan(Exp_time)]/5*Nt_Sgro
RampInput_Exp = Sgro2015Figure3excel["Ramp Input (nM cAMP)"]
RampInput_Exp = RampInput_Exp[~np.isnan(RampInput_Exp)]

Sgro_time = np.arange(0,16*Nt_Sgro,dt)
RampInput_Sgro= np.interp(Sgro_time, Exp_time,RampInput_Exp)

# test update function with one step caculations
# A,R,r=test_agent.update(1,dt)

# initializations
A_trace_orig=[A0]; R_trace_orig=[R0]; r_trace=[]

for i in range(len(Sgro_time)-1):
    A_now=A_trace_orig[i]
    R_now=R_trace_orig[i]
    signal_now=RampInput_Sgro[i]
    
    A_next,R_next,r_now=Sgro_agent.update(dt,signal_now)
    A_trace_orig.append(A_next)
    R_trace_orig.append(R_next)
    r_trace.append(r_now)
    
# Traces
A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
A_trace_plot=(A_trace_orig + Nh_Sgro_offset)/Nh_Sgro;
t_plot_Sgro = np.array(t)/Nt
# plot stimulus and traces
label_font_size = 10
trace_width = 3.0
#
fig5 = plt.figure(figsize=(6, 4))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Sgro,RampInput_Sgro, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
#ax1.set_title('cAMP stim:'+str(cAMP) + \
#              ', step stim:'+str(TStepOn)+'~'+str(TStepOff) + \
#              ', ramp stim:'+str(TRampOn)+'~'+str(TRampOff) )
ax1.set_title('cAMP stim:'+str(cAMP))
# ax1.set_xlabel('Time',fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot_Sgro,A_trace_plot, color='g',linewidth=trace_width)
ax2.set_ylabel('Activator',fontsize=label_font_size)
ax2.yaxis.label.set_color('g')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax2.set_xlabel('Time',fontsize=label_font_size)
plt.show()
#%%
#np.savetxt('Single_cell_ramp_Sgro_OUT_200218.txt', (t_plot_Sgro, A_trace_plot))
# read out the saved traces
[t_plot_Sgro, A_trace_plot] = np.loadtxt('Single_cell_ramp_Sgro_OUT_200218.txt', dtype=float)

#%% Goldbeter 1987
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

p0=0.8; a0=3; b0=0.9; g0=0

# Fig 2, 4 variable model, autonomous oscillation
Goldbeter3_agent=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
# initializations
p_trace=[p0]; b_trace=[b0]; g_trace=[g0]

# stimulation settings
cAMP=1
T_tot = 16; T_tot = T_tot*Nt_Goldbeter
dt=0.001; t=list(np.arange(0,T_tot,dt))

Exp_time = Sgro2015Figure3excel["Ramp Input (min Time)"]# np.linspace(0,20,num=1000) # t step 0.02
Exp_time = Exp_time[~np.isnan(Exp_time)]/5*Nt_Goldbeter
RampInput_Exp = Sgro2015Figure3excel["Ramp Input (nM cAMP)"]
RampInput_Exp = RampInput_Exp[~np.isnan(RampInput_Exp)]

Goldbeter_time = np.arange(0,16*Nt_Goldbeter,dt)
RampInput_Goldbeter= np.interp(Goldbeter_time, Exp_time,RampInput_Exp)

for i in range(len(t)-1):
    p_now=p_trace[i]
    signal_now = RampInput_Goldbeter[i]
    p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0, signal_now)
    p_trace.append(p_next)
    b_trace.append(b_next)
    g_trace.append(g_next)
        
   
# Convert into np array
b_trace = np.array(b_trace)/Nh_Goldbeter; 
p_trace = np.array(p_trace); 
t_plot_Goldbeter = np.array(t)/Nt_Goldbeter

# plot stimulus and traces
label_font_size = 10
trace_width = 3.0

fig5 = plt.figure(figsize=(6, 4))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Goldbeter,RampInput_Goldbeter, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
#ax1.set_title('cAMP stim:'+str(cAMP) + \
#              ', step stim:'+str(TStepOn)+'~'+str(TStepOff) + \
#              ', ramp stim:'+str(TRampOn)+'~'+str(TRampOff) )

# ax1.set_xlabel('Time',fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot_Goldbeter, b_trace, color='g',linewidth=trace_width)
# line2 = ax2.plot(t_plot_Goldbeter, g_trace, color='b',linewidth=trace_width)
ax2.set_ylabel('Activator',fontsize=label_font_size)
ax2.yaxis.label.set_color('g')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax2.set_xlabel('Time',fontsize=label_font_size)

plt.show()
#%% Maeda and Loomis 2004
from MaedaLoomis2004_agent_and_pop_FUN import MaedaLoomis2004_agent
# Maeda & Loomis 2004 parameters
k1=2.0; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=1.0; k8=1.3; k9=0.3; k10=0.8
k11=0.7; k12=4.9; k13=23; k14=4.5
MaedaAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}
ACA0=0.1; PKA0=0.1; ERK20=0.1; RegA0=0.1; cAMPi0=0.01; 
cAMPe0=0.1; CAR10=0.1
state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]
MaedaLoomis_agent=MaedaLoomis2004_agent([1,1],state0,MaedaAgentParam)

# initializations
ACA_trace=[ACA0]; PKA_trace=[PKA0]; ERK2_trace=[ERK20]
RegA_trace=[RegA0]; cAMPi_trace=[cAMPi0]; # cAMPe_trace=[cAMPe0]
CAR1_trace=[CAR10]

# stimulation settings
cAMP=1
T_tot = 16 * Nt_Maeda
dt=0.001; t=list(np.arange(0,T_tot,dt))

Exp_time = Sgro2015Figure3excel["Ramp Input (min Time)"]# np.linspace(0,20,num=1000) # t step 0.02
Exp_time = Exp_time[~np.isnan(Exp_time)]/5*Nt_Maeda
RampInput_Exp = Sgro2015Figure3excel["Ramp Input (nM cAMP)"]
RampInput_Exp = RampInput_Exp[~np.isnan(RampInput_Exp)]

Maeda_time = np.arange(0,16*Nt_Maeda,dt)
RampInput_Maeda= np.interp(Maeda_time, Exp_time,RampInput_Exp)

for i in range(len(t)-1):
    signal_now = RampInput_Maeda[i]
    ACA_now=ACA_trace[i]
    PKA_now=PKA_trace[i]
    ERK2_now=ERK2_trace[i]
    RegA_now=RegA_trace[i]
    cAMPi_now=cAMPi_trace[i]
    cAMPe_now=cAMPi_trace[i]
    CAR1_now=CAR1_trace[i]
    
    ACA_next,PKA_next,ERK2_next,RegA_next,\
    cAMPi_next,cAMPe_next,CAR1_next=MaedaLoomis_agent.update(dt,signal_now)
    
    ACA_trace.append(ACA_next)
    PKA_trace.append(PKA_next)
    ERK2_trace.append(ERK2_next)
    RegA_trace.append(RegA_next)
    cAMPi_trace.append(cAMPi_next)
    # cAMPe_trace.append(cAMPe_next)
    CAR1_trace.append(CAR1_next)
    

ERK2_trace = np.array(ERK2_trace) # vectorize p_trace
cAMPi_trace = np.array(cAMPi_trace)/Nh_Maeda
t_plot_Maeda = np.array(t)/Nt_Maeda

# plot stimulus and traces
label_font_size = 10
trace_width = 3.0

fig5 = plt.figure(figsize=(5, 4))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Maeda,RampInput_Maeda, linewidth=trace_width)
ax1.set_title('cAMP stim:'+str(cAMP))
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot_Maeda,cAMPi_trace, color='g',linewidth=trace_width)
ax2.set_ylabel(r'$cAMP_{cyto}$',fontsize=label_font_size)
ax2.yaxis.label.set_color('g')

ax3 = fig5.add_subplot(grid[1:, 0], sharex=ax2, frameon=False)
line2=ax3.plot(t_plot_Maeda,ERK2_trace, color='b', linewidth=trace_width) # right axis
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")
ax3.set_ylabel('ERK2' ,fontsize=label_font_size)
ax3.yaxis.label.set_color('b')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax3.set_xlabel('Time',fontsize=label_font_size)
plt.show()

#%% Kamino 2017
from Kamino2017_agent_and_pop_FUN import Kamino2017_agent 

tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
AgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
x0=0.01; y0=0.06; z0=0.005
Kamino_agent=Kamino2017_agent([x0,y0,z0],AgentParam)

x_trace=[x0]; y_trace=[y0]

# stimulation settings
cAMP=1
T_tot = 16 * Nt_Kamino
dt=0.001; t=list(np.arange(0,T_tot,dt))

Exp_time = Sgro2015Figure3excel["Ramp Input (min Time)"]# np.linspace(0,20,num=1000) # t step 0.02
Exp_time = Exp_time[~np.isnan(Exp_time)]/5*Nt_Kamino
RampInput_Exp = Sgro2015Figure3excel["Ramp Input (nM cAMP)"]
RampInput_Exp = RampInput_Exp[~np.isnan(RampInput_Exp)]

Kamino_time = np.arange(0,16*Nt_Kamino,dt)
RampInput_Kamino= np.interp(Kamino_time, Exp_time,RampInput_Exp)


for i in range(len(t)-1):
    signal_now = RampInput_Kamino[i]
    x_now=x_trace[i]
    y_now=y_trace[i]
    x_next,y_next,z_next= Kamino_agent.update(dt, signal_now)
    x_trace.append(x_next)
    y_trace.append(y_next)
        
   
# Convert into np array
x_trace = np.array(x_trace) # vectorize p_trace
y_trace = (np.array(y_trace)-Nh_Kamino_offset)/Nh_Kamino
t_plot_Kamino = np.array(t)/Nt_Kamino

# plot stimulus and traces
label_font_size = 10
trace_width = 3.0

fig5 = plt.figure(figsize=(5, 4))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Kamino,RampInput_Kamino, linewidth=trace_width)
ax1.set_title('cAMP stim:'+str(cAMP))
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot_Kamino,y_trace, color='g',linewidth=trace_width)
ax2.set_ylabel(r'$cAMP_{cyto}$',fontsize=label_font_size)
ax2.yaxis.label.set_color('g')

ax3 = fig5.add_subplot(grid[1:, 0], sharex=ax2, frameon=False)
line2=ax3.plot(t_plot_Kamino,x_trace, color='b', linewidth=trace_width) # right axis
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")
ax3.set_ylabel('x, adaptive component' ,fontsize=label_font_size)
ax3.yaxis.label.set_color('b')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax3.set_xlabel('Time',fontsize=label_font_size)
plt.show()


#%% Save all outputs in npz file
np.savez('single_cell_StepRamp_200320.npz', 
         t_plot_Goldbeter = t_plot_Goldbeter , b_trace=b_trace,
         t_plot_Maeda = t_plot_Maeda , cAMPi_trace = cAMPi_trace,
         t_plot_Sgro=t_plot_Sgro, A_trace_plot=A_trace_plot,
         t_plot_Kamino=t_plot_Kamino, y_trace=y_trace)
#%% load saved npz output file
npzfile = np.load('single_cell_StepRamp_200320.npz')
t_plot_Goldbeter =  npzfile['t_plot_Goldbeter'] ; b_trace = npzfile['b_trace']
t_plot_Maeda=npzfile['t_plot_Maeda'] ; cAMPi_trace=npzfile['cAMPi_trace']
t_plot_Sgro=npzfile['t_plot_Sgro']; A_trace_plot=npzfile['A_trace_plot']
t_plot_Kamino= npzfile['t_plot_Kamino']; y_trace=npzfile['y_trace']        


#%% Plot all model results- separate subplots
#title_font_size = 22
#label_font_size = 22
#tick_font_size = 16
#legend_font_size = 12
#trace_width = 3
abcd_font_size = 18
label_font_size=14
title_font_size = 14
sublabel_font_size = 12
trace_width=2
tick_font_size=12
mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']   

fig3 = plt.figure(figsize=(6,36))
grid = plt.GridSpec(11, 1, wspace=0.8, hspace=0.05)

ax0 = fig3.add_subplot(grid[0, 0])
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
ax0.text(-0.1, 1.4, 'A', ha='center',va='center',
         transform = ax0.transAxes, color = 'b', fontsize=abcd_font_size)

ax1= fig3.add_subplot(grid[1:2, 0])
ax1.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 1 FRET Trace"],
                               color='k', linewidth = trace_width)
ax1.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 2 FRET Trace"],
                               color='grey', linewidth = trace_width)
ax1.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 3 FRET Trace"],
                               color='lightgrey', linewidth = trace_width)
ax1.axvspan(10, 30, alpha=0.2, color='b'); ax1.axvspan(50, 70, alpha=0.2, color='b')
ax1.set_ylabel('FRET,\n A.U.',fontsize=sublabel_font_size)
ax1.set_xlabel('Time (min)',fontsize=sublabel_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_ylim([-0.1,0.6]); ax1.set_xlim([0,80])

ax2= fig3.add_subplot(grid[3:4, 0])
ax2.plot(t_plot_Goldbeter, b_trace, color=mycolors[0],linewidth=trace_width)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Martiel 1987',color = mycolors[0],  fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax2.axvspan(2, 6, alpha=0.3, color='g'); ax2.axvspan(10, 14, alpha=0.3, color='g')
ax2.set_ylim([-0.2,1.2]); ax2.set_xlim([0,16])
ax2.text(-0.1, 1.12, 'B', ha='center',va='center',
         transform = ax2.transAxes, color = 'g', fontsize=abcd_font_size)

ax3= fig3.add_subplot(grid[5:6, 0])
ax3.plot(t_plot_Maeda, cAMPi_trace, color=mycolors[1],linewidth=trace_width)
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Maeda 2004', color = mycolors[1],fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax3.axvspan(2,6, alpha=0.3, color='g'); ax3.axvspan(10,14, alpha=0.3, color='g')
ax3.set_ylim([-0.2,1.2]);  ax3.set_xlim([0,16])
ax3.text(-0.1, 1.12, 'C', ha='center',va='center',
         transform = ax3.transAxes, color = 'g', fontsize=abcd_font_size)

ax4= fig3.add_subplot(grid[7:8, 0])
ax4.plot(t_plot_Sgro, A_trace_plot, color=mycolors[5],linewidth=trace_width)
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('Sgro 2015', color = mycolors[5],fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax4.axvspan(2,6, alpha=0.3, color='g'); ax4.axvspan(10,14, alpha=0.3, color='g')
ax4.set_ylim([-0.4,1.2]);  ax4.set_xlim([0,16])
ax4.text(-0.1, 1.12, 'D', ha='center',va='center',
         transform = ax4.transAxes, color = 'g', fontsize=abcd_font_size)

ax5= fig3.add_subplot(grid[9:10, 0])
ax5.plot(t_plot_Kamino, y_trace, color=mycolors[7],linewidth=trace_width)
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5.set_title('Kamino 2017', color = mycolors[7],fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax5.axvspan(2,6, alpha=0.3, color='g'); ax5.axvspan(10,14, alpha=0.3, color='g')
ax5.set_ylim([-0.4,1.2]);  ax5.set_xlim([0,16])
ax5.text(-0.1, 1.12, 'E', ha='center',va='center',
         transform = ax5.transAxes, color = 'g', fontsize=abcd_font_size)

fig3.text(0.01, 0.5, r'$cAMP_{i}$', fontsize=label_font_size,va='center', rotation='vertical')
fig3.text(0.5, 0.12, 'Time, A.U.', fontsize=label_font_size, ha='center')
plt.show()


#%% Plot all model results- all one plot
title_font_size = 22
label_font_size = 22
tick_font_size = 16
legend_font_size = 12
trace_width = 3
mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']   

fig3 = plt.figure(figsize=(20,12))
grid = plt.GridSpec(3, 2, wspace=0.4, hspace=0.5)

ax0 = fig3.add_subplot(grid[0, 0])
ax0.plot(Sgro2015Figure3excel["Ramp Input (min Time)"],Sgro2015Figure3excel["Ramp Input (nM cAMP)"],
                              color='k', linewidth = trace_width)
ax0.set_ylabel(r'$cAMP_{ext}$(nM)',fontsize=label_font_size)
ax0.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax1= fig3.add_subplot(grid[1:, 0])
ax1.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 1 FRET Trace"],
                               color='k', linewidth = trace_width)
ax1.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 2 FRET Trace"],
                               color='dimgrey', linewidth = trace_width)
ax1.plot(Sgro2015Figure3excel["Cell Trace Time (min)"],Sgro2015Figure3excel["Cell 3 FRET Trace"],
                               color='darkgrey', linewidth = trace_width)
ax1.axvspan(10, 30, alpha=0.2, color='b'); ax1.axvspan(50, 70, alpha=0.2, color='b')
ax1.set_ylabel('FRET Signal, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Time (min)',fontsize=label_font_size)

ax2= fig3.add_subplot(grid[1:, 1])
ax2.plot(t_plot_Goldbeter, b_trace, color=mycolors[0],
          label='Goldbeter 1987', linewidth=trace_width)
ax2.plot(t_plot_Maeda, cAMPi_trace, color=mycolors[1],
         label = 'Maeda&Loomis 2004', linewidth=trace_width)
ax2.plot(t_plot_Sgro, A_trace_plot, color=mycolors[5],
         label = 'Sgro&Mehta 2015',linewidth=trace_width)
ax2.plot(t_plot_Kamino, y_trace, color=mycolors[7],
         label = 'Kamino 2017', linewidth=trace_width)
ax2.axvspan(2.5, 7.5, alpha=0.3, color='g'); ax2.axvspan(12.5, 17.5, alpha=0.3, color='g')

ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Model Results', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax2.set_ylabel(r'$cAMP_{cyto}$, A.U.',fontsize=label_font_size)
ax2.set_xlabel('Time, A.U.',fontsize=label_font_size)
leg = ax2.legend( frameon=False,loc='right center',ncol=1,prop={'size': 17})

plt.show()