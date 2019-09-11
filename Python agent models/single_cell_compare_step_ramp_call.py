# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:23:10 2019

@author: Chuqiao
"""

import numpy as np
import random
import math
import matplotlib.pyplot as plt

#%% 
Nt_Gregor = 6 
Nt_Sgro = 27
Nt_Goldbeter = 6.94
Nt_Maeda = 3.57
Nt_Maeda = 5.22
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

# stimulation settings
cAMP=1
T_tot = 22.5
dt=0.005 ; t_tot=T_tot*Nt_Sgro; t=list(np.arange(0,t_tot,dt))
signal_trace_Sgro=np.zeros(len(t))

# define step stimulus
TStepOn = 2.5; tStepOn = math.floor(TStepOn/dt*Nt); 
TStepOff = 7.5; tStepOff = math.floor(TStepOff/dt*Nt); 
signal_trace_Sgro[tStepOn:tStepOff] = cAMP

# define exponential ramp stimulus
TRampOn=10; tRampOn=math.floor(TRampOn/dt*Nt)
TRampOff=20; tRampOff=math.floor(TRampOff/dt*Nt)
tRamp=np.linspace(tRampOn, tRampOff, num=tRampOff-tRampOn+1)

a=math.log(cAMP+1)
temp_t=np.linspace(0,1,len(tRamp))
signal_trace_Sgro[tRampOn:tRampOff+1] = np.exp(a*temp_t)-1


# test update function with one step caculations
# A,R,r=test_agent.update(1,dt)

# initializations
A_trace_orig=[A0]; R_trace_orig=[R0]; r_trace=[]

for i in range(len(t)-1):
    A_now=A_trace_orig[i]
    R_now=R_trace_orig[i]
    signal_now=signal_trace_Sgro[i]
    
    A_next,R_next,r_now=Sgro_agent.update(dt,signal_now)
    A_trace_orig.append(A_next)
    R_trace_orig.append(R_next)
    r_trace.append(r_now)
    
# Traces
A_trace_offset=1.5
A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
A_trace_plot=(A_trace_orig+A_trace_offset)/Na;
t_plot_Sgro = np.array(t)/Nt
# plot stimulus and traces
label_font_size = 10
trace_width = 3.0

fig5 = plt.figure(figsize=(6, 4))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Sgro,signal_trace_Sgro, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
ax1.set_title('cAMP stim:'+str(cAMP) + \
              ', step stim:'+str(TStepOn)+'~'+str(TStepOff) + \
              ', ramp stim:'+str(TRampOn)+'~'+str(TRampOff) )

# ax1.set_xlabel('Time',fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot_Sgro,A_trace_plot, color='g',linewidth=trace_width)
ax2.set_ylabel('Activator',fontsize=label_font_size)
ax2.yaxis.label.set_color('g')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax2.set_xlabel('Time',fontsize=label_font_size)

plt.show()


#%% Goldbeter 1987
from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_agent_3var

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

# stimulation settings
cAMP=1
T_tot = 225
dt=0.001; t=list(np.arange(0,T_tot,dt))
signal_trace_Gold=np.zeros(len(t))

# define step stimulus
TStepOn = 25; tStepOn = math.floor(TStepOn/dt); 
TStepOff = 75; tStepOff = math.floor(TStepOff/dt); 
signal_trace_Gold[tStepOn:tStepOff] = cAMP

# define exponential ramp stimulus
TRampOn=100; tRampOn=math.floor(TRampOn/dt)
TRampOff=200; tRampOff=math.floor(TRampOff/dt)
tRamp=np.linspace(tRampOn, tRampOff, num=tRampOff-tRampOn+1)

a=math.log(cAMP+1)
temp_t=np.linspace(0,1,len(tRamp))
signal_trace_Gold[tRampOn:tRampOff+1] = np.exp(a*temp_t)-1


for i in range(len(t)-1):
    p_now=p_trace[i]
    p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0, signal_trace_Gold[i])
    p_trace.append(p_next)
    b_trace.append(b_next)
    g_trace.append(g_next)
        
   
# Convert into np array
b_trace = np.array(b_trace); # b_trace = b_trace/np.amax(b_trace)
p_trace = np.array(p_trace); # p_trace = p_trace/np.amax(p_trace)
t_plot_Goldbeter = np.array(t)

# plot stimulus and traces
label_font_size = 10
trace_width = 3.0

fig5 = plt.figure(figsize=(6, 4))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Goldbeter,signal_trace_Gold, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
ax1.set_title('cAMP stim:'+str(cAMP) + \
              ', step stim:'+str(TStepOn)+'~'+str(TStepOff) + \
              ', ramp stim:'+str(TRampOn)+'~'+str(TRampOff) )

# ax1.set_xlabel('Time',fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot_Goldbeter, b_trace, color='g',linewidth=trace_width)
# line2 = ax2.plot(t_plot_Goldbeter, g_trace, color='b',linewidth=trace_width)
ax2.set_ylabel('Activator',fontsize=label_font_size)
ax2.yaxis.label.set_color('g')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax2.set_xlabel('Time',fontsize=label_font_size)

plt.show()

''''
label_font_size=25; trace_width=3; tick_font_size=18
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
'''
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
T_tot = 225
dt=0.001; t=list(np.arange(0,T_tot,dt))
signal_trace_Maeda=np.zeros(len(t))

# define step stimulus
TStepOn = 25; tStepOn = math.floor(TStepOn/dt); 
TStepOff = 75; tStepOff = math.floor(TStepOff/dt); 
signal_trace_Maeda[tStepOn:tStepOff] = cAMP

# define exponential ramp stimulus
TRampOn=100; tRampOn=math.floor(TRampOn/dt)
TRampOff=200; tRampOff=math.floor(TRampOff/dt)
tRamp=np.linspace(tRampOn, tRampOff, num=tRampOff-tRampOn+1)

a=math.log(cAMP+1)
temp_t=np.linspace(0,1,len(tRamp))
signal_trace_Maeda[tRampOn:tRampOff+1] = np.exp(a*temp_t)-1


for i in range(len(t)-1):
    ACA_now=ACA_trace[i]
    PKA_now=PKA_trace[i]
    ERK2_now=ERK2_trace[i]
    RegA_now=RegA_trace[i]
    cAMPi_now=cAMPi_trace[i]
    cAMPe_now=cAMPi_trace[i]
    CAR1_now=CAR1_trace[i]
    
    ACA_next,PKA_next,ERK2_next,RegA_next,\
    cAMPi_next,cAMPe_next,CAR1_next=MaedaLoomis_agent.update(dt,signal_trace_Maeda[i])
    
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

# plot stimulus and traces
label_font_size = 10
trace_width = 3.0

fig5 = plt.figure(figsize=(5, 4))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Laub,signal_trace_Maeda, linewidth=trace_width)
ax1.set_title('cAMP stim:'+str(cAMP) + \
              ', step stim:'+str(TStepOn)+'~'+str(TStepOff) + \
              ', ramp stim:'+str(TRampOn)+'~'+str(TRampOff) )
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)

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

#%% Kamino 2017
from Kamino2017_agent_and_pop_FUN import Kamino2017_agent 

tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
AgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
x0=0.01; y0=0.05; z0=0.005
Kamino_agent=Kamino2017_agent([x0,y0,z0],AgentParam)

x_trace=[x0]; y_trace=[y0]

# stimulation settings
cAMP=1
T_tot = 225
dt=0.001; t=list(np.arange(0,T_tot,dt))
signal_trace_Kamino=np.zeros(len(t))

# define step stimulus
TStepOn = 25; tStepOn = math.floor(TStepOn/dt); 
TStepOff = 75; tStepOff = math.floor(TStepOff/dt); 
signal_trace_Kamino[tStepOn:tStepOff] = cAMP

# define exponential ramp stimulus
TRampOn=100; tRampOn=math.floor(TRampOn/dt)
TRampOff=200; tRampOff=math.floor(TRampOff/dt)
tRamp=np.linspace(tRampOn, tRampOff, num=tRampOff-tRampOn+1)

a=math.log(cAMP+1)
temp_t=np.linspace(0,1,len(tRamp))
signal_trace_Kamino[tRampOn:tRampOff+1] = np.exp(a*temp_t)-1


for i in range(len(t)-1):
    x_now=x_trace[i]
    y_now=y_trace[i]
    x_next,y_next,z_next= Kamino_agent.update(dt, signal_trace_Kamino[i])
    x_trace.append(x_next)
    y_trace.append(y_next)
        
   
# Convert into np array
x_trace = np.array(x_trace) # vectorize p_trace
y_trace = np.array(y_trace)
t_plot_Kamino = np.array(t)

# plot stimulus and traces
label_font_size = 10
trace_width = 3.0

fig5 = plt.figure(figsize=(5, 4))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot_Kamino,signal_trace_Kamino, linewidth=trace_width)
ax1.set_title('cAMP stim:'+str(cAMP) + \
              ', step stim:'+str(TStepOn)+'~'+str(TStepOff) + \
              ', ramp stim:'+str(TRampOn)+'~'+str(TRampOff) )
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

#%% PLot all the stuff
# Signal normalization
cAMPi_trace = cAMPi_trace/(np.amax(cAMPi_trace)) # Laub Loomis 1998
b_trace = b_trace/(np.amax(b_trace)) # Goldbeter 1987
# gregor_campCyto_trace=(gregor_campCyto_trace-np.amin(gregor_campCyto_trace))/np.amax(gregor_campCyto_trace-np.amin(gregor_campCyto_trace)) # Gregor2010
y_trace= (y_trace-np.amin(y_trace))/np.max((y_trace-np.amin(y_trace)))  # Kamino2017

A_trace_plot = A_trace_plot/(np.amax(A_trace_plot))
A_t_norm_param=5
A_trace_plot = A_trace_plot
t_plot_Sgro = t_plot_Sgro*A_t_norm_param # Sgro 2015

#%% Plottitle_font_size = 20 
label_font_size = 8
tick_font_size = 8
legend_font_size = 8
trace_width = 3
title_font_size = 12

#fig3 = plt.figure(figsize=(16, 14))
#grid = plt.GridSpec(2, 2, wspace=0.15, hspace=0.25)
fig3 = plt.figure(figsize=(3,13))
grid = plt.GridSpec(9, 1, wspace=0.05, hspace=1)

ax1= fig3.add_subplot(grid[0, 0])
ax1.plot(t_plot_Kamino,signal_trace_Kamino, linewidth=trace_width)
ax1.set_title('cAMP stimulation trace', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
# ax1.set_xlabel( 'Time, A.U.' ,fontsize=label_font_size)

ax2= fig3.add_subplot(grid[1:3, 0])
ax2.plot(t_plot_Goldbeter, b_trace, color='g',linewidth=trace_width)
ax2.set_ylabel(r'$cAMP_{cyto}$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Goldbeter 1987', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax3= fig3.add_subplot(grid[3:5, 0])
ax3.plot(t_plot_Laub, cAMPi_trace, color='g',linewidth=trace_width)
ax3.set_ylabel(r'$cAMP_{cyto}$',fontsize=label_font_size)
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Laub & Loomis 1998', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax4= fig3.add_subplot(grid[5:7, 0])
ax4.plot(t_plot_Sgro, A_trace_plot, color='g',linewidth=trace_width)
ax4.set_ylabel(r'$cAMP_{cyto}$',fontsize=label_font_size)
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('Sgro 2015', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax5= fig3.add_subplot(grid[7:9, 0])
ax5.plot(t_plot_Kamino, y_trace, color='g',linewidth=trace_width)
ax5.set_ylabel(r'$cAMP_{cyto}$',fontsize=label_font_size)
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5.set_title('Kamino 2017', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

plt.show()