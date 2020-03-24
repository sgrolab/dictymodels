# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:23:10 2019

@author: Chuqiao
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
#%% Adaptive spiking
# Gregor 2010
from Gregor2010_agent_and_pop_FUN import  Gregor2010_agent

Amax=20;  Abas=0.4 # uM
w=2*math.pi/6 # min-1
Vc=1.1e-9 # ml
St=1.33 # cm2
Sc=1.3e-6 # cm2
K=0.0004 # uM, 400 pM
c_sec= 3.6 # min-1
c_excite=1.01 # min-1

GregorAgentParam={'Amax':Amax,'Abas':Abas,'w':w,'Vc':Vc,'St':St,'Sc':Sc,'K':K,\
            'c_sec':c_sec,'c_excite':c_excite}
# Initializations
campCyto0 = 0.4
sinthetai0 = (campCyto0*2-Amax-Abas)/(-Amax+Abas)
thetai0 = np.arcsin(sinthetai0)
campExt0 = 0

Gregor_agent=Gregor2010_agent([1,1],[campCyto0, thetai0, campExt0],GregorAgentParam)

dt=0.001; t_tot=6*Nt_Gregor; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
#constant_signal=1e-6
constant_signal=1e-6#1e-6 # 
stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
signal_trace=np.zeros(len(t))
signal_trace[stim_time_step:] = constant_signal

eta=0.002 # noise
# initializations
gregor_thetai_trace=[thetai0]; gregor_campCyto_trace=[campCyto0]; gregor_r_trace=[0]

for i in range(len(t)-1):
    thetai_now=gregor_thetai_trace[i]
    campCyto_now=gregor_campCyto_trace[i]
    
    thetai_next, campCyto_next, r_now = Gregor_agent.update(dt,eta, signal_trace[i])
    gregor_thetai_trace.append( thetai_next)
    gregor_campCyto_trace.append(campCyto_next)
    gregor_r_trace.append(r_now)
    
#Traces
gregor_thetai_trace= np.array(gregor_thetai_trace) 
gregor_campCyto_trace= np.array(gregor_campCyto_trace) /Nh_Gregor
t_plot_Gregor = np.array(t)/Nt_Gregor
#t_plot_Gregor = np.array(t)*6/(Nt_Gregor)

## Use find peaks to get time and height narmalization factor
#PkPos, Props = find_peaks(gregor_campCyto_trace,  prominence=0.95)
#plt.plot(gregor_campCyto_trace)
#plt.plot(PkPos,gregor_campCyto_trace[PkPos], "x")
#Wdths = peak_widths(gregor_campCyto_trace, PkPos, rel_height=0.95)
#plt.hlines(*Wdths[1:], color="C2")
#Nt_Gregor = round(Wdths[0][0] * dt,2)
#Nh_Gregor = Props['prominences'][0]
#plt.title('For Gregor, time normalization parameter is: '+ str(Nt_Gregor)+
#      ' \n and height normalization parameter is: '+str(round(Nh_Gregor,2)))
#plt.show()

#  Plot adaptive spike
label_font_size = 25
trace_width = 6.0 # tim scale normalization factor based on intrinsic period

fig5 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Gregor,signal_trace, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
# ax1.set_xlabel('Time',fontsize=label_font_size)
ax1.set_title('cAMP stim level '+str(constant_signal))

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot_Gregor,gregor_campCyto_trace, color='g',linewidth=trace_width)
ax2.set_ylabel(r'$cAMP_{cyto}$',fontsize=label_font_size)
# ax2.yaxis.label.set_color('g')
ax2.set_xlabel('Time',fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot_Gregor,gregor_campCyto_trace, color='g',linewidth=trace_width)
ax2.set_ylabel(r'$cAMP_{cyto}$',fontsize=label_font_size)
# ax2.yaxis.label.set_color('g')
ax2.set_xlabel('Time',fontsize=label_font_size)

ax3 = fig5.add_subplot(grid[1:, 0], sharex=ax2, frameon=False)
line2=ax3.plot(t_plot_Gregor,gregor_thetai_trace, color='b', linewidth=trace_width) # right axis
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")
ax3.set_ylabel(r'$\theta$' ,fontsize=label_font_size)
ax3.yaxis.label.set_color('b')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))

plt.show()

#%% Sgro 2015
from Sgro2015_agent_and_pop_FUN import Sgro2015_agent

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}

A0=-1.5; R0=-0.5
Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)
dt=0.005 ; t_tot=6*Nt_Sgro; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
constant_signal_arr = np.array([1]) #np.logspace(-3, 2, num=5)
# constant_signal_arr = np.array([10000])
A_plot_all = np.zeros((len(constant_signal_arr),len(t)))

for j in range(len(constant_signal_arr)):
    constant_signal= constant_signal_arr[j]
    stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
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
    A_trace_plot=(A_trace_orig + Nh_Sgro_offset)/Nh_Sgro;
    A_plot_all[j,:] = A_trace_plot
    
## Check find_peaks
#PkPos, Props = find_peaks(A_trace_plot,  prominence=0.95)
#plt.plot(A_trace_plot)
#plt.plot(PkPos, A_trace_plot[PkPos], "x")
#Wdths = peak_widths(A_trace_plot, PkPos, rel_height=0.95)
#plt.hlines(*Wdths[1:], color="C2")
#Nt_Sgro = round(Wdths[0][0] * dt,2)
#Nh_Sgro = Props['prominences'][0]
#plt.title('For Sgro, time normalization parameter is: '+ str(Nt_Sgro)+
#      ' \n and height normalization parameter is: '+str(round(Nh_Sgro,2)))
#plt.show()

colors = plt.cm.summer(np.linspace(0,1,len(constant_signal_arr)+1))
t_plot_Sgro = np.array(t)/(Nt_Sgro)
fig3 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
ax1= fig3.add_subplot(grid[0, 0])

for j in range(len(constant_signal_arr))   :
    this_A_trace = A_plot_all[j,:]
    signal = constant_signal_arr[j]
    ax1.plot(t_plot_Sgro,this_A_trace,linewidth=2.5,
             color=colors[j], label = 'Input='+str(signal))
ax1.axvline(1,color='grey',linewidth='2.5')
ax1.set_ylabel( 'cAMPi trace, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax1.set_title('Sgro 2015')
leg = ax1.legend();

#%% Goldbeter 1987
from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_agent_3var
##Fig 5 parameters
#k1 = 0.036     # per min
#k2 = 0.666    # per min
#L1 = 10; L2 = 0.005 
#c = 10;           # 0.15 ~ 50
#lamda=0.01; theta=0.01
#e=  0.108 # compared to 1
#q=4000
#sig= 0.57 # compared to 0.6
#v=12; k= 4 # k prime in the paper
#ki=0.958 # compared to 1.7 
#kt=0.9
#kc=3.58 # compared to 5.4
#h=5

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

dt=0.001; t_tot=6*Nt_Goldbeter; t=list(np.arange(0,t_tot,dt))

stim_time_step=int(round(1/6*t_tot/dt)) # at this time step step input is applied
#constant_signal=1
# constant_signal=10000

constant_signal_arr = np.array([1]) #np.logspace(-3, 1, num=5)
# constant_signal_arr = np.array([10000])
b_plot_all = np.zeros((len(constant_signal_arr),len(t)))

for j in range(len(constant_signal_arr)):
    constant_signal= constant_signal_arr[j]
    stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
    signal_trace=np.zeros(len(t))
    signal_trace[stim_time_step:] = constant_signal
    # initializations
    p_trace=[p0]; b_trace=[b0]; g_trace=[g0]

    for i in range(len(t)-1):
        p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,signal_trace[i])
        p_trace.append(p_next)
        b_trace.append(b_next)
        g_trace.append(g_next)
           
    # Convert into np array
    b_trace = np.array(b_trace)/Nh_Goldbeter; #b_trace = b_trace/np.amax(b_trace)
    p_trace = np.array(p_trace); #p_trace = p_trace/np.amax(p_trace)
    b_plot_all[j]=b_trace

## Check find_peaks
#PkPos, Props = find_peaks(b_trace,  prominence=0.95)
#plt.plot(b_trace)
#plt.plot(PkPos, b_trace[PkPos], "x")
#Wdths = peak_widths(b_trace, PkPos, rel_height=0.95)
#plt.hlines(*Wdths[1:], color="C2")
#Nt_Goldbeter = round(Wdths[0][0] * dt,2)
#Nh_Goldbeter = Props['prominences'][0]
#plt.title('For Goldbeter, time normalization parameter is: '+ str(Nt_Goldbeter)+
#      '\n and height normalization parameter is: '+str(round(Nh_Goldbeter,2)))
#plt.show()

colors = plt.cm.summer(np.linspace(0,1,len(constant_signal_arr)+1))
t_plot_Goldbeter = np.array(t)/(Nt_Goldbeter)
fig3 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
ax1= fig3.add_subplot(grid[0, 0])

for j in range(len(constant_signal_arr))   :
    this_b_trace = b_plot_all[j,:]
    signal = constant_signal_arr[j]
    ax1.plot(t_plot_Goldbeter,this_b_trace,linewidth=2.5,
             color=colors[j], label = 'Input='+str(signal))
ax1.axvline(1,color='grey',linewidth='2.5')
ax1.set_ylabel( 'cAMPi trace, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax1.set_title('Martiel 1987')
ax1.set_xlim([0,6])
leg = ax1.legend();
#%%
# Check find_peaks
#peaks, properties = find_peaks(b_trace,   prominence=0.95)
#plt.plot(b_trace)
#plt.plot(peaks, b_trace[peaks], "x")
#Wdths = peak_widths(b_trace, peaks, rel_height=0.95)
#plt.hlines(*Wdths[1:], color="C2")
#plt.show()
#Nt_Goldbeter = round(Wdths[0][0] * dt,2)
t_plot_Goldbeter = np.array(t)*6/(Nt_Goldbeter)

label_font_size=25; trace_width=5; tick_font_size=18

fig5 = plt.figure(figsize=(10, 12))
grid = plt.GridSpec(3, 1, wspace=0.2, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Goldbeter,signal_trace, linewidth=trace_width)

# ax1.set_ylim([-0.1,1.1*np.amax(signal_trace)])
ax1.set_ylabel( r'$cAMP_{ext}$'+'\n input, A.U.' ,fontsize=label_font_size)
ax1.get_xaxis().set_visible(False)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
# ax1.set_xlabel('Time',fontsize=label_font_size)
ax1.set_title('Oscillation after step cAMP input, different parameter regime')

ax2= fig5.add_subplot(grid[1:, 0])
ax2.plot(t_plot_Goldbeter, b_trace,linewidth=trace_width, label='Goldbeter 1987, oscillation parameters')

ax2.set_ylim([-0.2,1.3])
ax2.set_ylabel(r'$cAMP_{cyto}$, A.U.',fontsize=label_font_size)
ax2.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)

leg = ax2.legend();
ax2.legend( frameon=False,loc='upper center',ncol=3,prop={'size': 15})
# ax2.legend((line1, line2), ('Gregor 2010','Sgro & Mehta 2015'))

plt.show()
#%% Maeda Loomis 2004
from MaedaLoomis2004_agent_and_pop_FUN import MaedaLoomis2004_agent

## Laub loomis parameters
#k1=1.4; k2=0.9; k3=2.5; k4=1.5; k5=0.6
#k6=0.8; k7=2.0; k8=1.3; k9=0.7; k10=1.0
#k11=0.3; k12=3.1; k13=1.8; k14=1.5

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
Maeda_agent=MaedaLoomis2004_agent([1,1],state0,MaedaAgentParam)

# Simulate time traces
"""
# defining signals
cAMP=1
tTrig=int(round(1/dt*Nt)) # Time step at which triggering happens
signal_trace=[0]*len(t)
signal_trace[tTrig:]=[cAMP]*len(signal_trace[tTrig:])

# test update function with one step caculations
# A,R,r=test_agent.update(1,dt)
"""
# initializations
ACA_trace=[ACA0]; PKA_trace=[PKA0]; ERK2_trace=[ERK20]
RegA_trace=[RegA0]; cAMPi_trace=[cAMPi0]; # cAMPe_trace=[cAMPe0]
CAR1_trace=[CAR10]

dt=0.0001; t_tot=6*Nt_Maeda; t=list(np.arange(0,t_tot,dt))

stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
##constant_signal=1
#constant_signal=1000
#signal_trace=np.zeros(len(t))
#signal_trace[stim_time_step:] = constant_signal

constant_signal_arr = np.array([1]) #np.logspace(-3, 1, num=5)
cAMPi_plot_all = np.zeros((len(constant_signal_arr),len(t)))

for j in range(len(constant_signal_arr)):
    constant_signal= constant_signal_arr[j]
    stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
    signal_trace=np.zeros(len(t))
    signal_trace[stim_time_step:] = constant_signal
    # initializations
    ACA_trace=[ACA0]; PKA_trace=[PKA0]; ERK2_trace=[ERK20]
    RegA_trace=[RegA0]; cAMPi_trace=[cAMPi0]; # cAMPe_trace=[cAMPe0]
    CAR1_trace=[CAR10]
    for i in range(len(t)-1):
        ACA_now=ACA_trace[i]
        PKA_now=PKA_trace[i]
        ERK2_now=ERK2_trace[i]
        RegA_now=RegA_trace[i]
        cAMPi_now=cAMPi_trace[i]
        cAMPe_now=cAMPi_trace[i]
        CAR1_now=CAR1_trace[i]
        
        ACA_next,PKA_next,ERK2_next,RegA_next,\
        cAMPi_next,cAMPe_next,CAR1_next=Maeda_agent.update(dt,signal_trace[i])
        
        ACA_trace.append(ACA_next)
        PKA_trace.append(PKA_next)
        ERK2_trace.append(ERK2_next)
        RegA_trace.append(RegA_next)
        cAMPi_trace.append(cAMPi_next)
        # cAMPe_trace.append(cAMPe_next)
        CAR1_trace.append(CAR1_next)
        
    
    ERK2_trace = np.array(ERK2_trace) # vectorize p_trace
    cAMPi_trace = np.array(cAMPi_trace)/Nh_Maeda
    cAMPi_plot_all[j,:] = cAMPi_trace
    
## Check find_peaks
#PkPos, Props = find_peaks(cAMPi_trace,  prominence=0.95)
#plt.plot(cAMPi_trace)
#plt.plot(PkPos, cAMPi_trace[PkPos], "x")
#Wdths = peak_widths(cAMPi_trace, PkPos, rel_height=0.95)
#plt.hlines(*Wdths[1:], color="C2")
#Nt_Maeda = round(Wdths[0][0] * dt,2)
#Nh_Maeda = Props['prominences'][0]
#plt.title('For Maeda, time normalization parameter is: '+ str(Nt_Maeda)+
#      '\n and height normalization parameter is: '+str(round(Nh_Maeda,2)))
#plt.show()

colors = plt.cm.summer(np.linspace(0,1,len(constant_signal_arr)+1))
t_plot_Maeda = np.array(t)/(Nt_Maeda)
fig3 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
ax1= fig3.add_subplot(grid[0, 0])

for j in range(len(constant_signal_arr))   :
    this_cAMPi_trace = cAMPi_plot_all[j,:]
    signal = constant_signal_arr[j]
    ax1.plot(t_plot_Maeda,this_cAMPi_trace,linewidth=2.5,
             color=colors[j], label = 'Input='+str(signal))
ax1.axvline(1,color='grey',linewidth='2.5')
ax1.set_ylabel( 'cAMPi trace, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax1.set_title('Maeda 2004')
ax1.set_xlim([0,6])
leg = ax1.legend();   
    
    #%%
## Check find_peaks
#peaks, properties = find_peaks(cAMPi_trace,   prominence=1)
#plt.plot(cAMPi_trace)
#plt.plot(peaks, cAMPi_trace[peaks], "x")
#Wdths = peak_widths(cAMPi_trace, peaks, rel_height=0.95)
#plt.hlines(*Wdths[1:], color="C2")
#plt.show()
#Nt_Maeda = round(Wdths[0][0] * dt,2)

t_plot_Maeda = np.array(t)*6/(Nt_Maeda)
# Figure single cell adaptive spike 
label_font_size = 25

trace_width = 6.0

fig5 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Maeda,signal_trace, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
# ax1.set_xlabel('Time',fontsize=label_font_size)
ax1.set_title('cAMP stim:'+str(constant_signal))

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

#%% Adaptive spike for individual cells
from Kamino2017_agent_and_pop_FUN import Kamino2017_agent 

tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
AgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
x0=0.01; y0=0.06; z0=0.005
Kamino_agent=Kamino2017_agent([x0,y0,z0],AgentParam)

x_trace=[x0]; y_trace=[y0]

dt=0.0001; t_tot=6*Nt_Kamino; t=list(np.arange(0,t_tot,dt))

#stim_time_step1=int(round(1/6*t_tot/dt)) # at this time step input is applied
#stim_time_step2=int(round(0.6*t_tot/dt)) 
##constant_signal1=1; constant_signal2=1
#constant_signal1=10000; constant_signal2=10000
#signal_trace=np.zeros(len(t))
#signal_trace[stim_time_step1:stim_time_step2] = constant_signal1
#signal_trace[stim_time_step2:] = constant_signal2

stim_time_step=int(round(1/6*t_tot/dt)) # at this time step step input is applied
constant_signal_arr = np.array([1]) #np.logspace(-5, 0, num=6)#############
y_plot_all = np.zeros((len(constant_signal_arr),len(t)))
for j in range(len(constant_signal_arr)):
    constant_signal= constant_signal_arr[j]
    stim_time_step=int(round(1/6*t_tot/dt)) # at this time step input is applied
    signal_trace=np.zeros(len(t))
    signal_trace[stim_time_step:] = constant_signal
    # initializations
    x_trace=[x0]; y_trace=[y0]
    for i in range(len(t)-1):
        x_now=x_trace[i]
        y_now=y_trace[i]
        x_next,y_next,z_next= Kamino_agent.update( dt, signal_trace[i])
        x_trace.append(x_next)
        y_trace.append(y_next)
               
    # Convert into np array
    x_trace = np.array(x_trace) 
    y_trace = (np.array(y_trace)- Nh_Kamino_offset) / Nh_Kamino
    y_plot_all[j,:] = y_trace

## Check find_peaks and find out the Nt and Nh
#PkPos, Props = find_peaks(y_trace,  prominence=0.2)
#plt.plot(y_trace)
#plt.plot(PkPos, y_trace[PkPos], "x")
#Wdths = peak_widths(y_trace, PkPos, rel_height=0.95)
#plt.hlines(*Wdths[1:], color="C2")
#Nt_Kamino = round(Wdths[0][0] * dt,2)
#Nh_Kamino= Props['prominences'][0]
#plt.title('For Kamino, time normalization parameter is: '+ str(Nt_Kamino)+
#      '\n and height normalization parameter is: '+str(round(Nh_Kamino,2)))
#plt.show()

colors = plt.cm.summer(np.linspace(0,1,len(constant_signal_arr)+1))
t_plot_Kamino = np.array(t)/(Nt_Kamino)
fig3 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
ax1= fig3.add_subplot(grid[0, 0])

for j in range(len(constant_signal_arr))   :
    this_y_trace = y_plot_all[j,:]
    signal = constant_signal_arr[j]
    ax1.plot(t_plot_Kamino,this_y_trace,linewidth=2.5,
             color=colors[j], label = 'Input='+str(signal))
ax1.axvline(1,color='grey',linewidth='2.5')
ax1.set_ylabel( 'cAMPi trace, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax1.set_title('Kamino 2017')
ax1.set_xlim([0,6])
leg = ax1.legend();
#%%
## Check find_peaks
#peaks, properties = find_peaks(y_trace,   prominence=0.2)
#plt.plot(y_trace)
#plt.plot(peaks, y_trace[peaks], "x")
#Wdths = peak_widths(y_trace, peaks, rel_height=0.95)
#plt.hlines(*Wdths[1:], color="C2")
#plt.show()
#Nt_Kamino = round(Wdths[0][0] * dt,2)
#print(Nt_Kamino)

t_plot_Kamino = np.array(t)*6/Nt_Kamino

#Plot adaptive spike
label_font_size = 15
trace_width = 6

fig5 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(2, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Kamino ,signal_trace, linewidth=trace_width)
# ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
# ax1.set_xlabel('Time',fontsize=label_font_size)
# ax1.set_title('cAMP stim:'+str(constant_signal))

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot_Kamino ,y_trace, color='g',linewidth=trace_width)
# ax2.set_ylabel(r'$cAMP_{cyto}$, y',fontsize=label_font_size)
# ax2.yaxis.label.set_color('g')

#ax3 = fig5.add_subplot(grid[1:, 0], sharex=ax2, frameon=False)
#line2=ax3.plot(t_plot_Kamino ,x_trace, color='b', linewidth=trace_width) # right axis
#ax3.yaxis.tick_right()
#ax3.yaxis.set_label_position("right")
#ax3.set_ylabel('Adaptation variable, x' ,fontsize=label_font_size)
#ax3.yaxis.label.set_color('b')
## ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
#ax3.set_xlabel('Time',fontsize=label_font_size)
plt.show()

## Signal normalization
#cAMPi_trace = cAMPi_trace/(np.amax(cAMPi_trace)) # Maeda Loomis 1998
#b_trace = b_trace/(np.amax(b_trace)) # Goldbeter 1987
#gregor_campCyto_trace=(gregor_campCyto_trace-np.amin(gregor_campCyto_trace))/np.amax(gregor_campCyto_trace-np.amin(gregor_campCyto_trace)) # # for sustained oscillation plot
##gregor_campCyto_trace=gregor_campCyto_trace-np.amin(gregor_campCyto_trace) # for adaptive spike plot
#y_trace= (y_trace-np.amin(y_trace))/np.max((y_trace-np.amin(y_trace)))  # Kamino2017
#A_trace_plot = A_trace_plot/(np.amax(A_trace_plot))
#%% Save all outputs in npz file
np.savez('single_cell_adaptive_200319.npz', 
         t_plot_Goldbeter = t_plot_Goldbeter , b_trace=b_trace,
         t_plot_Maeda=t_plot_Maeda, cAMPi_trace=cAMPi_trace,
         t_plot_Gregor=t_plot_Gregor, gregor_campCyto_trace=gregor_campCyto_trace,
         t_plot_Sgro=t_plot_Sgro, A_trace_plot=A_trace_plot,
         t_plot_Kamino=t_plot_Kamino, y_trace=y_trace)
         
#%% experimental data
my_dir = r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/'
Sgro2015Figure1excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheetname='Figure1')
#%% load saved npz output file
# npzfile = np.load('single_cell_oscillation_200319.npz')
npzfile = np.load('single_cell_adaptive_200319.npz')
t_plot_Goldbeter =  npzfile['t_plot_Goldbeter'] ; b_trace = npzfile['b_trace']
t_plot_Maeda=npzfile['t_plot_Maeda'] ; cAMPi_trace=npzfile['cAMPi_trace']
t_plot_Gregor=npzfile['t_plot_Gregor']; gregor_campCyto_trace=npzfile['gregor_campCyto_trace']
t_plot_Sgro=npzfile['t_plot_Sgro']; A_trace_plot=npzfile['A_trace_plot']
t_plot_Kamino= npzfile['t_plot_Kamino']; y_trace=npzfile['y_trace']
#%% Plot single cell adaptive spikes
abcd_font_size = 28
label_font_size=24
title_font_size = 26
sublabel_font_size = 22
trace_width=3
tick_font_size=20
mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']          
fig5 = plt.figure(figsize=(9.5, 8))
grid = plt.GridSpec(3, 1, wspace=1.2, hspace=1.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 1 FRET Trace (1nM)"],
                              linewidth=trace_width,color='k')
ax1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 2 FRET Trace (1nM)"],
                               linewidth=trace_width,color='dimgrey')
ax1.plot(Sgro2015Figure1excel["Time (min)"],Sgro2015Figure1excel["Cell 3 FRET Trace (1nM)"],
                               linewidth=trace_width,color='darkgrey')
#ax1.set_ylabel(r'FRET Signal, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Time (min)',fontsize=label_font_size)
ax1.axvline(x=5, ls='--', linewidth=trace_width, color=mycolors[6]) #dashed line at 5 (cAMP onset)
ax1.set_ylim([-0.1,0.7]); ax1.set_xlim([0, 30])
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Experiment',color = 'k', fontsize = title_font_size)
ax1.text(0.7,0.8,' 10uM cAMP', horizontalalignment='center',verticalalignment='center',
         transform = ax1.transAxes, color = 'k', fontsize=label_font_size)
ax1.text(-0.13, 1.5, 'A',
         horizontalalignment='center',verticalalignment='center',
         transform = ax1.transAxes, color = 'b', fontsize=abcd_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
#ax2.axvline(x = 5 , ymin=-0.2, ymax = 1.2, ls = '--', 
#            linewidth=trace_width, color=mycolors[6])
ax2.plot(t_plot_Goldbeter, b_trace,linewidth=trace_width, color=mycolors[0],
         label='Martiel 1987')
# ax3 = ax2.twinx()
ax3= fig5.add_subplot(grid[1:, 0])
ax3.plot(t_plot_Maeda, cAMPi_trace,linewidth=trace_width,color=mycolors[1],
         label='Maeda 2004')
#ax3.set_ylabel(r'MAeda 2004 $cAMP_{i}$, A.U.', color=mycolors[1],
#               fontsize=label_font_size-3)
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
# ax3.set_ylim([-30,260])

ax2.plot(t_plot_Gregor,gregor_campCyto_trace,linewidth=trace_width, 
         color=mycolors[2],label='Gregor 2010')
ax2.plot(t_plot_Sgro, A_trace_plot,linewidth=trace_width,
         color=mycolors[5], label='Sgro 2015')
ax2.plot(t_plot_Kamino, y_trace,linewidth=trace_width, 
         ls = '--',color=mycolors[7], label='Kamino 2017')
ax2.set_title('Simulation',color = 'k', fontsize = title_font_size)
ax2.text(-0.13, 1.06, 'B',
         horizontalalignment='center',verticalalignment='center',
         transform = ax2.transAxes, color = 'g', fontsize=abcd_font_size)

# ax2.set_ylim([-0.2,1.85])# oscillations
ax2.set_ylim([-0.2,1.2]) # adaptive spikes
ax2.set_xlim([0,6])
ax2.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
#########set y tick to 2 digits
fig5.text(0.03, 0.82, 'FRET, A.U.', ha='center', va='center',rotation='vertical',
          fontsize=label_font_size)
fig5.text(0.03, 0.35,r'$cAMP_{i}$, A.U.', ha='center', va='center',rotation='vertical',
          fontsize=label_font_size)
#leg = ax2.legend()
#ax2.legend( frameon=False,loc='upper center',ncol=1,prop={'size': 17})
plt.show()