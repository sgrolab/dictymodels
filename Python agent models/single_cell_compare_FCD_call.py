# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:23:10 2019

@author: Chuqiao
"""
# Fold change detection
import numpy as np
import random
import math
import matplotlib.pyplot as plt

from scipy import signal
from scipy.signal import find_peaks

#%% Kamino 2017
from Kamino_agent_single_cell import Kamino2017_agent 

z0First_space_Kamino = np.array([0.1, 1,3,10]) # first period extracellular cAMP level # np.array([100]) 
FC_space_Kamino= np.logspace(0.1, 2, num=8) # np.logspace(0.5, 0.5, num=1) 

# Initialize 
PkPrm_Kamino = np.zeros((len(z0First_space_Kamino), len(FC_space_Kamino))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus

tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
AgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
x0=0.01; y0=0.05; z0=0.005
Kamino_agent=Kamino2017_agent([x0,y0,z0],AgentParam)


dt=0.001; t_tot=100; t=list(np.arange(0,t_tot,dt))

for j in range(len(z0First_space_Kamino)):
    signal_trace=z0First_space_Kamino[j]*np.ones(len(t))
    for k in range(len(FC_space_Kamino)):
        
        stim_time_step=int(round(0.5*t_tot/dt)) # at this time second step input is applied
        signal_trace[stim_time_step:] = FC_space_Kamino[k]*z0First_space_Kamino[j]
        x_trace=[x0]; y_trace=[y0]
        
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
        
#        # check traces
#        plt.figure()
#        # plt.plot(t_plot_Kamino,x_trace,t_plot_Kamino,y_trace,t_plot_Kamino,signal_trace)
#        plt.plot(t_plot_Kamino,y_trace)
#        plt.xlabel('Time')
#        plt.ylabel('x,y,z')
#        plt.title('cAMP from  '+str(z0First_space[j])+' to '+str(FC_space[k]*z0First_space[j]))
#        # plt.gca().legend(('x','y','z'))
#        plt.show()
        
        y_trace_second=y_trace[stim_time_step:] # the second part of trace, second spike
        PkPos, PkProperties = find_peaks(y_trace_second, prominence=(0,10))
        # Check find_peaks
        # plt.plot(z_trace_later)
        # plt.plot(peaks, z_trace_later[peaks], "x")
        if PkPos.size: # if there is a second spike
            PkPrm_Kamino[j,k]=PkProperties["prominences"][0]
        else:
            PkPrm_Kamino[j,k]=0 # if there is no second spike
 
           
# plot FC vs. second response amplitude
label_font_size = 9
trace_width = 2

fig3 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)

ax1= fig3.add_subplot(grid[0, 0])

for i in range(len(z0First_space_Kamino)):
    ax1.plot(FC_space_Kamino,PkPrm_Kamino[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Kamino[i]))

ax1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
ax1.set_xlabel('Fold change in extracellular cAMP',fontsize=label_font_size)
ax1.set_xscale('log')
ax1.set_title('Second spike prominence')
leg = ax1.legend();

#%% Sgro 2015, w/o noise (sigma=0)
from Sgro2015_agent import Sgro2015_agent

z0First_space_Sgro =np.array([0.5,1,2]) # np.array([0.125,0.25,0.5,1,2,4]) # first period extracellular cAMP level # 
FC_space_Sgro=  np.logspace(0.5, 1.6, num=8) # np.logspace(0.5, 0.5, num=1) 
run_time_space =[0]# np.arange(0,0,1) # run 10 times and plot the mean of each case
# Initialize 

PkPrm_Sgro_woNoise = np.zeros((len(z0First_space_Sgro), len(FC_space_Sgro), len(run_time_space))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
Nt=27; # normalization factor of t
Na=3.5;  # normalization factor of A
A0=-1.5; R0=-0.5
Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)

dt=0.001 ; t_tot=20*Nt; t=list(np.arange(0,t_tot,dt))

# j_test=[2]; k_test=[3]

for j in range(len(z0First_space_Sgro)):
    signal_trace=z0First_space_Sgro[j]*np.ones(len(t))
    for k in range(len(FC_space_Sgro)):
        
        stim_time_step=int(round(0.5*t_tot/dt)) # at this time second step input is applied
        signal_trace[stim_time_step:] = FC_space_Sgro[k]*z0First_space_Sgro[j]
        for test in run_time_space:
            A_trace_orig=[A0]; R_trace_orig=[R0]
            for i in range(len(t)-1):
                A_now=A_trace_orig[i]
                R_now=R_trace_orig[i]
                signal_now=signal_trace[i]
        
                A_next,R_next,r_now=Sgro_agent.update(signal_now,dt)
                A_trace_orig.append(A_next)
                R_trace_orig.append(R_next)
       
            A_trace_offset=1.5
            A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
            A_trace_plot=(A_trace_orig+A_trace_offset)/Na;
            t_plot_Sgro = np.array(t)/Nt
            
    #        # check traces
    #        fig3 = plt.figure(figsize=(6, 6))
    #        grid = plt.GridSpec(2, 1, wspace=0.3, hspace=0.2)
    #        ax1= fig3.add_subplot(grid[0, 0])
    #        ax1.plot(t_plot_Sgro,signal_trace)
    #        ax1.set_ylabel('extracellular cAMP')
    #        ax1.set_title('cAMP from '+str(z0First_space_Sgro[j])+' to FC '+ str(FC_space_Sgro[k]))
    #        ax2= fig3.add_subplot(grid[1, 0])
    #        ax2.plot(t_plot_Sgro,A_trace_plot)
    #        ax2.set_ylabel('Activator')
    #        ax2.set_ylabel('Time')
    #        plt.show()
    
            end_time_step=stim_time_step + int(round(2.5*Nt/dt))
            A_trace_second=A_trace_plot[stim_time_step:end_time_step] # the second part of trace, second spike
            PkPos, PkProperties = find_peaks(A_trace_second, prominence=(0.4,2))
            # Check find_peaks
            # plt.plot(z_trace_later)
            # plt.plot(peaks, z_trace_later[peaks], "x")
            if PkPos.size: # if there is a second spike
                PkPrm_Sgro_woNoise[j,k,test]=PkProperties["prominences"][0]
            else:
                PkPrm_Sgro_woNoise[j,k,test]=0 # if there is no second spike
 
PkPrm_Sgro_mean=np.mean(PkPrm_Sgro_woNoise,axis=2)
# plot FC vs. second response amplitude
fig3 = plt.figure(figsize=(16, 14))
grid = plt.GridSpec(2, 2, wspace=0.1, hspace=0.35)

title_font_size = 20 
label_font_size = 20
tick_font_size = 16
legend_font_size = 15
trace_width = 5
ax3= fig3.add_subplot(grid[1, 0])
for i in range(len(z0First_space_Sgro)):
    ax3.plot(FC_space_Sgro,PkPrm_Sgro_woNoise[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Sgro[i]))

ax3.set_ylim([-0.3,1])
ax3.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
ax3.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
# ax3.set_xscale('log')
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Sgro 2015', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax3.legend();
ax3.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})

#%% Sgro 2015
from Sgro2015_agent import Sgro2015_agent

z0First_space_Sgro_noise =np.array([0.5,1,2]) # np.array([0.125,0.25,0.5,1,2,4]) # first period extracellular cAMP level # 
FC_space_Sgro_noise=  np.logspace(0.5, 1.6, num=8) # np.logspace(0.5, 0.5, num=1) 
run_time_space =np.arange(0,10,1) # run 10 times and plot the mean of each case
# Initialize 

PkPrm_Sgro_noise = np.zeros((len(z0First_space_Sgro_noise), len(FC_space_Sgro_noise), len(run_time_space))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
Nt=27; # normalization factor of t
Na=3.5;  # normalization factor of A
A0=-1.5; R0=-0.5
Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)

dt=0.001 ; t_tot=20*Nt; t=list(np.arange(0,t_tot,dt))

# j_test=[2]; k_test=[3]

for j in range(len(z0First_space_Sgro_noise)):
    signal_trace=z0First_space_Sgro_noise[j]*np.ones(len(t))
    for k in range(len(FC_space_Sgro_noise)):
        
        stim_time_step=int(round(0.5*t_tot/dt)) # at this time second step input is applied
        signal_trace[stim_time_step:] = FC_space_Sgro_noise[k]*z0First_space_Sgro_noise[j]
        for test in run_time_space:
            A_trace_orig=[A0]; R_trace_orig=[R0]
            for i in range(len(t)-1):
                A_now=A_trace_orig[i]
                R_now=R_trace_orig[i]
                signal_now=signal_trace[i]
        
                A_next,R_next,r_now=Sgro_agent.update(signal_now,dt)
                A_trace_orig.append(A_next)
                R_trace_orig.append(R_next)
       
            A_trace_offset=1.5
            A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
            A_trace_plot=(A_trace_orig+A_trace_offset)/Na;
            t_plot_Sgro_noise = np.array(t)/Nt
            
    #        # check traces
    #        fig3 = plt.figure(figsize=(6, 6))
    #        grid = plt.GridSpec(2, 1, wspace=0.3, hspace=0.2)
    #        ax1= fig3.add_subplot(grid[0, 0])
    #        ax1.plot(t_plot_Sgro,signal_trace)
    #        ax1.set_ylabel('extracellular cAMP')
    #        ax1.set_title('cAMP from '+str(z0First_space_Sgro[j])+' to FC '+ str(FC_space_Sgro[k]))
    #        ax2= fig3.add_subplot(grid[1, 0])
    #        ax2.plot(t_plot_Sgro,A_trace_plot)
    #        ax2.set_ylabel('Activator')
    #        ax2.set_ylabel('Time')
    #        plt.show()
    
            end_time_step=stim_time_step + int(round(2.5*Nt/dt))
            A_trace_second=A_trace_plot[stim_time_step:end_time_step] # the second part of trace, second spike
            PkPos, PkProperties = find_peaks(A_trace_second, prominence=(0.4,2))
            # Check find_peaks
            # plt.plot(z_trace_later)
            # plt.plot(peaks, z_trace_later[peaks], "x")
            if PkPos.size: # if there is a second spike
                PkPrm_Sgro_noise[j,k,test]=PkProperties["prominences"][0]
            else:
                PkPrm_Sgro_noise[j,k,test]=0 # if there is no second spike
 
PkPrm_Sgro_mean_noise=np.mean(PkPrm_Sgro_noise,axis=2)
# plot FC vs. second response amplitude
fig3 = plt.figure(figsize=(16, 14))
grid = plt.GridSpec(2, 2, wspace=0.1, hspace=0.35)

title_font_size = 20 
label_font_size = 20
tick_font_size = 16
legend_font_size = 15
trace_width = 5
ax3= fig3.add_subplot(grid[1, 0])
for i in range(len(z0First_space_Sgro)):
    ax3.plot(FC_space_Sgro_noise,PkPrm_Sgro_mean_noise[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Sgro_noise[i]))

ax3.set_ylim([-0.3,1])
ax3.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
ax3.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax3.set_xscale('log')
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Sgro 2015 (with noise)', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax3.legend();
ax3.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})
#%% Goldbeter 1987
from Goldbeter1987_agent import Goldbeter1987_agent_3var

z0First_space_Gold =np.array([0.1,0.2,0.4,0.8]) #np.array([0.5,10]) # first period extracellular cAMP level # 
FC_space_Gold=  np.logspace(0.5, 2, num=8) # np.logspace(0.5, 0.5, num=1) 

# Initialize 
PkPrm_Gold = np.zeros((len(z0First_space_Gold), len(FC_space_Gold))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus

k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10; L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01; theta=0.01
e= 0.108 # compared to 1
q=4000
sig=0.57 # compared to 0.6
v=12; k= 4 # k prime in the paper
ki=0.958 # compared to 1.7 
kt=0.9
kc=3.58 # compared to 5.4
h=5

Goldbeter3AgentParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}

p0=0.8; a0=3; b0=0.9; g0=0

dt=0.001; t_tot=100; t=list(np.arange(0,t_tot,dt))

# j_test=[3]; k_test=np.array([7])

for j in range(len(z0First_space_Gold)):
    signal_trace=z0First_space_Gold[j]*np.ones(len(t))
    for k in range(len(FC_space_Gold)):
        
        stim_time_step=int(round(0.5*t_tot/dt)) # at this time second step input is applied
        signal_trace[stim_time_step:] = FC_space_Gold[k]*z0First_space_Gold[j]
        
        p_trace=[p0]; b_trace=[b0]; g_trace=[g0]
        Goldbeter3_agent=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
        for i in range(len(t)-1):
            p_now=p_trace[i]
            p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,signal_trace[i])
            p_trace.append(p_next)
            b_trace.append(b_next)
            g_trace.append(g_next)
   
        b_trace = np.array(b_trace); b_trace = b_trace/np.amax(b_trace)
        t_plot_Gold = np.array(t)
        
#        # check traces
#        fig3 = plt.figure(figsize=(6, 6))
#        grid = plt.GridSpec(2, 1, wspace=0.3, hspace=0.2)
#        ax1= fig3.add_subplot(grid[0, 0])
#        ax1.plot(t_plot_Gold,signal_trace)
#        ax1.set_ylabel('extracellular cAMP')
#        ax1.set_title('cAMP from '+str(z0First_space_Gold[j])+' to FC '+ str(FC_space_Gold[k]))
#        ax2= fig3.add_subplot(grid[1, 0])
#        ax2.plot(t_plot_Gold,b_trace)
#        ax2.set_ylabel('beta, [cAMP]cyto')
#        ax2.set_ylabel('Time')
#        plt.show()

        end_time_step=len(signal_trace)
        b_trace_second=b_trace[stim_time_step:end_time_step] # the second part of trace, second spike
        PkPos, PkProperties = find_peaks(b_trace_second, prominence=(0,2))
        # Check find_peaks
        # plt.plot(z_trace_later)
        # plt.plot(peaks, z_trace_later[peaks], "x")
        if PkPos.size: # if there is a second spike
            PkPrm_Gold[j,k]=PkProperties["prominences"][0]
        else:
            PkPrm_Gold[j,k]=0 # if there is no second spike
                

# PkPrm_Gold_mean=np.mean(PkPrm_Gold,axis=2)
# plot FC vs. second response amplitude
label_font_size = 9
trace_width = 2

fig3 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)

ax1= fig3.add_subplot(grid[0, 0])

for i in range(len(z0First_space_Gold)):
    ax1.plot(FC_space_Gold,PkPrm_Gold[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Gold[i]))

ax1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
ax1.set_xscale('log')
ax1.set_xlabel('Fold change in extracellular cAMP',fontsize=label_font_size)
ax1.set_title('second spike prominence')
leg = ax1.legend();

        
   
#%% Laub and Loomis 1998
from LaubLoomis1998_agent import LaubLoomis1998_agent

z0First_space_Loomis =np.array([0.1,0.2,0.4, 0.8]) #np.array([0.5,10]) # first period extracellular cAMP level # 
FC_space_Loomis=  np.logspace(0.4, 1.6, num=8) # np.logspace(0.5, 0.5, num=1) 

# Initialize 
PkPrm_Loomis = np.zeros((len(z0First_space_Loomis), len(FC_space_Loomis))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus

k1=1.4; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=2.0; k8=1.3; k9=0.7; k10=1.0
k11=0.3; k12=3.1; k13=1.8; k14=1.5
LaubAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}
ACA0=0.1; PKA0=0.1; ERK20=0.1; RegA0=0.1; cAMPi0=0.01; 
cAMPe0=0.1; CAR10=0.1
state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]

dt=0.0005; t_tot=100; t=list(np.arange(0,t_tot,dt))

#  k_test=np.array([7])
j_test=[0,1,2,3,4,5];k_test=[0,1,2,3,4,5,6,7]

for j in range(len(z0First_space_Loomis)):  # j_test:#
    signal_trace=z0First_space_Loomis[j]*np.ones(len(t))
    for k in range(len(FC_space_Gold)):# k_test: #
        
        stim_time_step=int(round(0.5*t_tot/dt)) # at this time second step input is applied
        signal_trace[stim_time_step:] = FC_space_Loomis[k]*z0First_space_Loomis[j]
        
        ACA_trace=[ACA0]; PKA_trace=[PKA0]; ERK2_trace=[ERK20]; RegA_trace=[RegA0]; 
        cAMPi_trace=[cAMPi0]; cAMPe_trace=[cAMPe0]; CAR1_trace=[CAR10]
        
        LaubLoomis_agent=LaubLoomis1998_agent([1,1],state0,LaubAgentParam)
        
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
        t_plot_Loomis = np.array(t)
#        # check traces
#        fig3 = plt.figure(figsize=(3, 3))
#        grid = plt.GridSpec(2, 1, wspace=0.3, hspace=0.2)
#        ax1= fig3.add_subplot(grid[0, 0])
#        ax1.plot(t_plot_Loomis,signal_trace)
#        ax1.set_ylabel('extracellular cAMP')
#        ax1.set_title('cAMP from '+str(z0First_space_Loomis[j])+' to FC '+ str(FC_space_Loomis[k]))
#        ax2= fig3.add_subplot(grid[1, 0])
#        ax2.plot(t_plot_Loomis,cAMPi_trace)
#        ax2.set_ylim([0,25])
#        ax2.set_ylabel(' [cAMP]cyto')
#        ax2.set_ylabel('Time')
#        plt.show()

        end_time_step=len(signal_trace)
        cAMPi_trace_second=cAMPi_trace[stim_time_step:end_time_step] # the second part of trace, second spike
        PkPos, PkProperties = find_peaks(cAMPi_trace_second, prominence=(0,20))
        # Check find_peaks
        # plt.plot(z_trace_later)
        # plt.plot(peaks, z_trace_later[peaks], "x")
        if PkPos.size: # if there is a second spike
            PkPrm_Loomis[j,k]=PkProperties["prominences"][0]
        else:
            PkPrm_Loomis[j,k]=0 # if there is no second spike
                

# PkPrm_Gold_mean=np.mean(PkPrm_Gold,axis=2)
# plot FC vs. second response amplitude
label_font_size = 9
trace_width = 2

fig3 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)

ax1= fig3.add_subplot(grid[0, 0])

for i in range(len(z0First_space_Loomis)):
    ax1.plot(FC_space_Loomis,PkPrm_Loomis[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Loomis[i]))

ax1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
ax1.set_xscale('log')
ax1.set_xlabel('Fold change in extracellular cAMP',fontsize=label_font_size)
ax1.set_title('second spike prominence')
leg = ax1.legend();

#%% Plot all 4 models
title_font_size = 20 
label_font_size = 20
tick_font_size = 16
legend_font_size = 15
trace_width = 5

#fig3 = plt.figure(figsize=(16, 14))
#grid = plt.GridSpec(2, 2, wspace=0.15, hspace=0.25)
fig3 = plt.figure(figsize=(12,13))
grid = plt.GridSpec(2, 2, wspace=0.1, hspace=0.3)

ax1= fig3.add_subplot(grid[0, 0])
for i in range(len(z0First_space_Gold)):
    ax1.plot(FC_space_Gold,PkPrm_Gold[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Gold[i]))

ax1.set_ylim([-0.3,1])
#ax1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#ax1.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax1.set_xscale('log')
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Goldbeter 1987', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax1.legend();
ax1.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})

ax2= fig3.add_subplot(grid[0, 1])
for i in range(len(z0First_space_Loomis)):
    ax2.plot(FC_space_Loomis,PkPrm_Loomis[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Loomis[i]))

ax2.set_ylim([-1,3])
#ax2.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#ax2.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax2.set_xscale('log')
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Laub & Loomis 1998', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax2.legend();
ax2.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})


ax3= fig3.add_subplot(grid[1, 0])
for i in range(len(z0First_space_Sgro)):
    ax3.plot(FC_space_Sgro_noise,PkPrm_Sgro_mean_noise[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Sgro_noise[i]))

ax3.set_ylim([-0.3,0.9])
#ax3.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#ax3.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax3.set_xscale('log')
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Sgro 2015 (with noise)', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax3.legend();
ax3.legend( frameon=False,loc='bottom',ncol=2,prop={'size': legend_font_size})


ax4= fig3.add_subplot(grid[1, 1])
for i in range(len(z0First_space_Kamino)):
    ax4.plot(FC_space_Kamino,PkPrm_Kamino[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Kamino[i]))
ax4.set_ylim([-0.08,0.27])    
#ax4.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#ax4.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax4.set_xscale('log')
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('Kamino & Sawai 2017', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax4.legend();
ax4.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})

fig3.text(0.5, 0.04, r'$cAMP_{ext}$'+' fold change', ha='center', va='center',fontsize=label_font_size)
fig3.text(0.06, 0.5, 'Second spike prominence', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
#%% Plot all 4 models

title_font_size = 20 
label_font_size = 20
tick_font_size = 16
legend_font_size = 14
trace_width = 5

#fig3 = plt.figure(figsize=(16, 14))
#grid = plt.GridSpec(2, 2, wspace=0.15, hspace=0.25)
fig3 = plt.figure(figsize=(12,13))
grid = plt.GridSpec(2, 2, wspace=0.2, hspace=0.2)

ax1= fig3.add_subplot(grid[0, 0])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Gold)) )
for i in range(len(z0First_space_Gold)):
    ax1.plot(FC_space_Gold,PkPrm_Gold[i,:], color=colors[i],linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Gold[i]))


ax1.set_ylim([-0.3,1])
#ax1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#ax1.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax1.set_xscale('log')
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Goldbeter 1987', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax1.legend();
ax1.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})

ax2= fig3.add_subplot(grid[0, 1])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Loomis)) )
for i in range(len(z0First_space_Loomis)):
    ax2.plot(FC_space_Loomis,PkPrm_Loomis[i,:],color=colors[i], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Loomis[i]))

ax2.set_ylim([-1,3])
#ax2.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#ax2.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax2.set_xscale('log')
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Laub & Loomis 1998', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax2.legend();
ax2.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})


ax3= fig3.add_subplot(grid[1, 0])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Sgro)) )
for i in range(len(z0First_space_Sgro)):
    ax3.plot(FC_space_Sgro_noise,PkPrm_Sgro_woNoise[i,:],color=colors[i], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Sgro_noise[i]))

ax3.set_ylim([-0.3,0.9])
#ax3.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#ax3.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax3.set_xscale('log')
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Sgro 2015 (w/o noise)', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax3.legend();
ax3.legend( frameon=False,loc='bottom',ncol=2,prop={'size': legend_font_size})


ax4= fig3.add_subplot(grid[1, 1])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Kamino)) )
for i in range(len(z0First_space_Kamino)):
    ax4.plot(FC_space_Kamino,PkPrm_Kamino[i,:],color=colors[i], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Kamino[i]))
ax4.set_ylim([-0.08,0.27])    
#ax4.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#ax4.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax4.set_xscale('log')
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('Kamino & Sawai 2017', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax4.legend();
ax4.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})

fig3.text(0.5, 0.04, r'$cAMP_{ext}$'+' fold change', ha='center', va='center',fontsize=label_font_size)
fig3.text(0.06, 0.5, 'Second spike prominence', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
 