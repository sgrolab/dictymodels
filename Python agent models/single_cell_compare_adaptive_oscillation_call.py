# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:23:10 2019

@author: Chuqiao
"""

import numpy as np
import random
import math
import matplotlib.pyplot as plt

#%% Adaptive spiking
# Gregor 2010
from Gregor2010_agent import *

Amax=20;  Abas=0.4 # uM
w=2*pi/6 # min-1
Vc=1.1e-9 # ml
St=1.33 # cm2
Sc=1.3e-6 # cm2
K=0.0004 # uM, 400 pM
c_sec= 3.6 # min-1
c_excite=1.01 # min-1

GregorAgentParam={'Amax':Amax,'Abas':Abas,'w':w,'Vc':Vc,'St':St,'Sc':Sc,'K':K,\
            'c_sec':c_sec,'c_excite':c_excite}

campCyto0 = 0.4
sinthetai0 = (campCyto0*2-Amax-Abas)/(-Amax+Abas)
thetai0 = np.arcsin(sinthetai0)
campExt0 = 0

Gregor_agent=Gregor2010_agent([1,1],[campCyto0, thetai0, campExt0],GregorAgentParam)

dt=0.00005; t_tot=100; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
constant_signal=1e-6 # randomly set 2019/7/17
stim_time_step=int(round(0.2*t_tot/dt)) # at this time step input is applied
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
gregor_campCyto_trace= np.array(gregor_campCyto_trace) 
t_plot_Gregor = np.array(t)/(t_tot/25)

# Plot adaptive spike
label_font_size = 25
trace_width = 6.0

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
from Sgro2015_agent import *

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
Nt=27; # normalization factor of t
Na=3.5;  # normalization factor of A
A0=-1.5; R0=-0.5
Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)

dt=0.001 ; t_tot=5*Nt; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
constant_signal=1; # randomly set 2019/7/17
stim_time_step=int(round(0.2*t_tot/dt)) # at this time step input is applied
signal_trace=np.zeros(len(t))
signal_trace[stim_time_step:] = constant_signal


# test update function with one step caculations
# A,R,r=test_agent.update(1,dt)

# initializations
A_trace_orig=[A0]; R_trace_orig=[R0]; r_trace=[]

for i in range(len(t)-1):
    A_now=A_trace_orig[i]
    R_now=R_trace_orig[i]
    signal_now=signal_trace[i]
    
    A_next,R_next,r_now=Sgro_agent.update(signal_now,dt)
    A_trace_orig.append(A_next)
    R_trace_orig.append(R_next)
    r_trace.append(r_now)
    
# Traces
A_trace_offset=1.5
A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
A_trace_plot=(A_trace_orig+A_trace_offset)/Na;
t_plot_Sgro = np.array(t)/Nt

label_font_size = 10
trace_width = 6.0

fig5 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Sgro,signal_trace, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
ax1.set_title('cAMP stim:'+str(constant_signal))
# ax1.set_xlabel('Time',fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot_Sgro,A_trace_plot, color='g',linewidth=trace_width)
ax2.set_ylabel('Activator',fontsize=label_font_size)
ax2.yaxis.label.set_color('g')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax2.set_xlabel('Time',fontsize=label_font_size)

plt.show()

#%% Goldbeter 1987
from Goldbeter1987_agent import Goldbeter1987_agent_3var

k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10; L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01; theta=0.01
e=  0.108 # compared to 1
q=4000
sig= 0.57 # compared to 0.6
v=12; k= 4 # k prime in the paper
ki=0.958 # compared to 1.7 
kt=0.9
kc=3.58 # compared to 5.4
h=5

Goldbeter3AgentParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}

p0=0.8; a0=3; b0=0.9; g0=0

# Fig 2, 4 variable model, autonomous oscillation
Goldbeter3_agent=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
# initializations
p_trace=[p0]; b_trace=[b0]; g_trace=[g0]

dt=0.00005; t_tot=40; t=list(np.arange(0,t_tot,dt))

stim_time_step=int(round(0.5*t_tot/dt)) # at this time step step input is applied
constant_signal=1
signal_trace=np.zeros(len(t))
signal_trace[stim_time_step:] = constant_signal


for i in range(len(t)-1):
    p_now=p_trace[i]
    if i< stim_time_step:
        p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,'none')
    else:
        p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,constant_signal)
    p_trace.append(p_next)
    b_trace.append(b_next)
    g_trace.append(g_next)
        
   
# Convert into np array
b_trace = np.array(b_trace); b_trace = b_trace/np.amax(b_trace)
p_trace = np.array(p_trace); p_trace = p_trace/np.amax(p_trace)
t_plot_Goldbeter = np.array(t)

#label_font_size = 25; trace_width = 6.0
#fig5 = plt.figure(figsize=(6, 6))
#grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
#
#ax1= fig5.add_subplot(grid[0, 0])
#ax1.plot(t_plot_Goldbeter,signal_trace, linewidth=trace_width)
#ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
#ax1.set_title('cAMP stim:'+str(constant_signal))
## ax1.set_xlabel('Time',fontsize=label_font_size)
#
#ax2= fig5.add_subplot(grid[1:, 0])
#ax2.plot(t_plot_Goldbeter,b_trace, color='g',linewidth=trace_width, label='beta')
#ax2.plot(t_plot_Goldbeter,p_trace, color='b',linewidth=trace_width,label='pT')
#ax2.set_ylabel(r'$cAMP_{cyto},  \beta$',fontsize=label_font_size)
#ax2.yaxis.label.set_color('g')
#ax2.set_xlabel('Time',fontsize=label_font_size)
#leg = ax2.legend();
#ax2.legend( frameon=False,loc='right',ncol=1,prop={'size': 15})
#plt.show()


label_font_size=25; trace_width=5; tick_font_size=18

fig5 = plt.figure(figsize=(10, 12))
grid = plt.GridSpec(3, 1, wspace=0.2, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Goldbeter,signal_trace, linewidth=trace_width)

# ax1.set_ylim([-0.1,1.1*np.amax(signal_trace)])
ax1.set_xlim([0,t_tot])
ax1.set_ylabel( r'$cAMP_{ext}$'+'\n input, A.U.' ,fontsize=label_font_size)
ax1.get_xaxis().set_visible(False)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
# ax1.set_xlabel('Time',fontsize=label_font_size)
ax1.set_title('Oscillation after step cAMP input, different parameter regime')

ax2= fig5.add_subplot(grid[1:, 0])
ax2.plot(t_plot_Goldbeter, b_trace,linewidth=trace_width, label='Goldbeter 1987, oscillation parameters')

ax2.set_ylim([-0.2,1.3])
ax2.set_xlim([0,25])
ax2.set_ylabel(r'$cAMP_{cyto}$, A.U.',fontsize=label_font_size)
ax2.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)

leg = ax2.legend();
ax2.legend( frameon=False,loc='upper center',ncol=3,prop={'size': 15})
# ax2.legend((line1, line2), ('Gregor 2010','Sgro & Mehta 2015'))

plt.show()
#%% Laub and Loomis 1998
from LaubLoomis1998_agent import LaubLoomis1998_agent

k1=1.4; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=2.0; k8=1.3; k9=0.7; k10=1.0
k11=0.3; k12=3.1; k13=1.8; k14=1.5
LaubAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}
ACA0=0.1; PKA0=0.1; ERK20=0.1; RegA0=0.1; cAMPi0=0.01; 
cAMPe0=0.1; CAR10=0.1
state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]
LaubLoomis_agent=LaubLoomis1998_agent([1,1],state0,LaubAgentParam)

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

dt=0.00005; t_tot=25; t=list(np.arange(0,t_tot,dt))

stim_time_step=int(round(0.2*t_tot/dt)) # at this time step input is applied
constant_signal=1 # randomly set 2019/7/17
signal_trace=np.zeros(len(t))
signal_trace[stim_time_step:] = constant_signal

for i in range(len(t)-1):
    ACA_now=ACA_trace[i]
    PKA_now=PKA_trace[i]
    ERK2_now=ERK2_trace[i]
    RegA_now=RegA_trace[i]
    cAMPi_now=cAMPi_trace[i]
    cAMPe_now=cAMPi_trace[i]
    CAR1_now=CAR1_trace[i]
    
    ACA_next,PKA_next,ERK2_next,RegA_next,\
    cAMPi_next,cAMPe_next,CAR1_next=LaubLoomis_agent.update(dt,signal_trace[i])
    
    ACA_trace.append(ACA_next)
    PKA_trace.append(PKA_next)
    ERK2_trace.append(ERK2_next)
    RegA_trace.append(RegA_next)
    cAMPi_trace.append(cAMPi_next)
    # cAMPe_trace.append(cAMPe_next)
    CAR1_trace.append(CAR1_next)
    

ERK2_trace = np.array(ERK2_trace) # vectorize p_trace
cAMPi_trace = np.array(cAMPi_trace)
t_plot_Laub = np.array(t)


# Figure single cell adaptive spike 
label_font_size = 25
trace_width = 6.0

fig5 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Laub,signal_trace, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
# ax1.set_xlabel('Time',fontsize=label_font_size)
ax1.set_title('cAMP stim:'+str(constant_signal))

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot_Laub,cAMPi_trace, color='g',linewidth=trace_width)
ax2.set_ylabel(r'$cAMP_{cyto}$',fontsize=label_font_size)
ax2.yaxis.label.set_color('g')

ax3 = fig5.add_subplot(grid[1:, 0], sharex=ax2, frameon=False)
line2=ax3.plot(t_plot_Laub,ERK2_trace, color='b', linewidth=trace_width) # right axis
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")
ax3.set_ylabel('ERK2' ,fontsize=label_font_size)
ax3.yaxis.label.set_color('b')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax3.set_xlabel('Time',fontsize=label_font_size)
plt.show()

#%% Adaptive spike for individual cells
from Kamino_agent_single_cell import Kamino2017_agent 

tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
AgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
x0=0.01; y0=0.05; z0=0.005
Kamino_agent=Kamino2017_agent([x0,y0,z0],AgentParam)

x_trace=[x0]; y_trace=[y0]

dt=0.001; t_tot=100; t=list(np.arange(0,t_tot,dt))

#stim_time_step=int(round(0.2*t_tot/dt)) # at this time step input is applied
#constant_signal=1
#signal_trace=np.zeros(len(t))
#signal_trace[stim_time_step:] = constant_signal

stim_time_step1=int(round(0.2*t_tot/dt)) # at this time step input is applied
stim_time_step2=int(round(0.6*t_tot/dt)) 
constant_signal1=1; constant_signal2=3
signal_trace=np.zeros(len(t))
signal_trace[stim_time_step1:stim_time_step2] = constant_signal1
signal_trace[stim_time_step2:] = constant_signal2


for i in range(len(t)-1):
    x_now=x_trace[i]
    y_now=y_trace[i]
    x_next,y_next,z_next= Kamino_agent.update(1, dt, signal_trace[i])
    x_trace.append(x_next)
    y_trace.append(y_next)
        
   
# Convert into np array
x_trace = np.array(x_trace) # vectorize p_trace
y_trace = np.array(y_trace)
t_plot_Kamino = np.array(t)

#label_font_size = 25; trace_width = 6.0
#fig5 = plt.figure(figsize=(6, 6))
#grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
#
#ax1= fig5.add_subplot(grid[0, 0])
#ax1.plot(t_plot_Kamino,signal_trace, linewidth=trace_width)
#ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
#ax1.set_title('cAMP stim:'+str(constant_signal))
## ax1.set_xlabel('Time',fontsize=label_font_size)
#
#ax2= fig5.add_subplot(grid[1:, 0])
#ax2.plot(t_plot_Kamino,x_trace, color='g',linewidth=trace_width, label='inh. x')
#ax2.plot(t_plot_Kamino,y_trace, color='b',linewidth=trace_width,label='y')
#ax2.set_ylabel(r'$cAMP_{cyto}, y$',fontsize=label_font_size)
#ax2.yaxis.label.set_color('g')
#ax2.set_xlabel('Time',fontsize=label_font_size)
#leg = ax2.legend();
#ax2.legend( frameon=False,loc='right',ncol=1,prop={'size': 15})
#plt.show()

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
#%% Plot adaptive spike
# Signal normalization
cAMPi_trace = cAMPi_trace/(np.amax(cAMPi_trace)) # Laub Loomis 1998
b_trace = b_trace/(np.amax(b_trace)) # Goldbeter 1987
gregor_campCyto_trace=(gregor_campCyto_trace-np.amin(gregor_campCyto_trace))/np.amax(gregor_campCyto_trace-np.amin(gregor_campCyto_trace)) # Gregor2010
y_trace= (y_trace-np.amin(y_trace))/np.max((y_trace-np.amin(y_trace)))  # Kamino2017

A_trace_plot = A_trace_plot/(np.amax(A_trace_plot))
A_t_norm_param=5
A_trace_plot = A_trace_plot
t_plot_Sgro = t_plot_Sgro*A_t_norm_param # Sgro 2015

#%% Plot single cell adaptive spikes
t_tot=25
dt=0.1
t_plot_signal= np.arange(0,t_tot,dt)
stim_time_step=int(round(0.2*t_tot/dt)) # at this time step input is applied
constant_signal=1
signal_trace_plot=np.zeros(len(t_plot_signal))
signal_trace_plot[stim_time_step:] = constant_signal
 
label_font_size=25
trace_width=5
tick_font_size=18

fig5 = plt.figure(figsize=(10, 12))
grid = plt.GridSpec(3, 1, wspace=0.2, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_signal,signal_trace_plot, linewidth=trace_width)

ax1.set_ylim([-0.1,1.1*np.amax(signal_trace_plot)])
ax1.set_xlim([0,t_tot])
ax1.set_ylabel( r'$cAMP_{ext}$'+'\n input, A.U.' ,fontsize=label_font_size)
ax1.get_xaxis().set_visible(False)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
# ax1.set_xlabel('Time',fontsize=label_font_size)
ax1.set_title('Adaptive spike after step cAMP input')

ax2= fig5.add_subplot(grid[1:, 0])
ax2.plot(t_plot_Goldbeter, b_trace,linewidth=trace_width, label='Goldbeter 1987')
ax2.plot(t_plot_Laub, cAMPi_trace,linewidth=trace_width, label='Laub & Loomis 1998')
ax2.plot(t_plot_Gregor,gregor_campCyto_trace,linewidth=trace_width, label='Gregor 2010')
ax2.plot(t_plot_Sgro, A_trace_plot,linewidth=trace_width, label='Sgro & Mehta 2010')
ax2.plot(t_plot_Kamino, y_trace,linewidth=trace_width, label='Kamino & Sawai 2017')

ax2.set_ylim([-0.2,1.3])
ax2.set_xlim([0,25])
ax2.set_ylabel(r'$cAMP_{cyto}$, A.U.',fontsize=label_font_size)
ax2.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)

leg = ax2.legend();
ax2.legend( frameon=False,loc='upper center',ncol=3,prop={'size': 15})
# ax2.legend((line1, line2), ('Gregor 2010','Sgro & Mehta 2015'))

plt.show()