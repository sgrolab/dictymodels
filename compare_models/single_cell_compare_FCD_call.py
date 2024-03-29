# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:23:10 2019

@author: Chuqiao
"""
# Fold change detection
import pandas as pd
import numpy as np
import random
import math
import pandas as pd
import matplotlib.pyplot as plt
# set matplotlib default font
import matplotlib
font = {'family' : 'Roboto'}
matplotlib.rc('font', **font)

from scipy import signal
from scipy.signal import find_peaks
import scipy.io

# Normalization parameters
from NormParam import *
#%% Kamino 2017
from Kamino2017_agent_and_pop_FUN import Kamino2017_agent 

z0First_space_Kamino = np.array([1,2,4,8]) # first period extracellular cAMP level # np.array([100]) 
FC_space_Kamino= np.logspace(0.1, 1.6, num=8) # np.logspace(0.5, 0.5, num=1) 

# Initialize 
PkPrm_Kamino = np.zeros((len(z0First_space_Kamino), len(FC_space_Kamino))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus
PkPrm_Kamino_norm = np.zeros((len(z0First_space_Kamino), len(FC_space_Kamino)))

tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
AgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
x0=0.01; y0=0.06; z0=0.005
Kamino_agent=Kamino2017_agent([x0,y0,z0],AgentParam)


dt=0.001; t_tot=30 * Nt_Kamino; t=list(np.arange(0,t_tot,dt))

for j in range(len(z0First_space_Kamino)):
    
    signal_trace=z0First_space_Kamino[j]*np.zeros(len(t))
    
    for k in range(len(FC_space_Kamino)):
        stim_time_step1 = int(round(1/3*t_tot/dt))
        stim_time_step2=int(round(2/3*t_tot/dt)) # at this time second step input is applied
        
        signal_trace[stim_time_step1:stim_time_step2] = z0First_space_Kamino[j]
        signal_trace[stim_time_step2:] = FC_space_Kamino[k]*z0First_space_Kamino[j]
        x_trace=[x0]; y_trace=[y0]
        
        for i in range(len(t)-1):
            x_now=x_trace[i]
            y_now=y_trace[i]
            x_next,y_next,z_next= Kamino_agent.update( dt, signal_trace[i])
            x_trace.append(x_next)
            y_trace.append(y_next)
   
        # Convert into np array
        x_trace = np.array(x_trace) # vectorize p_trace
        y_trace = (np.array(y_trace)-Nh_Kamino_offset)/Nh_Kamino
        t_plot_Kamino = np.array(t)
        
#        # check traces
#        plt.figure(figsize=(4,2.5))
#        # plt.plot(t_plot_Kamino,x_trace,t_plot_Kamino,y_trace,t_plot_Kamino,signal_trace)
#        plt.plot(t_plot_Kamino/Nt_Kamino,y_trace)
#        plt.xlabel('Time')
#        plt.ylabel('x,y,z')
#        plt.title('cAMP from  '+str(z0First_space_Kamino[j])+' to '+str(FC_space_Kamino[k]*z0First_space_Kamino[j]))
#        # plt.gca().legend(('x','y','z'))
#        plt.show()
        
        
        # normalize 2nd peak height to the first peak height 
        y_trace_first=y_trace[stim_time_step1:stim_time_step2] 
        y_trace_second=y_trace[stim_time_step2:] # the second part of trace, second spike
        PkPos1, PkProperties1 = find_peaks(y_trace_first, prominence=(0,50/Nh_Kamino))
        PkPos2, PkProperties2 = find_peaks(y_trace_second, prominence=(0,50/Nh_Kamino))
        # Check find_peaks
        # plt.plot(z_trace_later)
        # plt.plot(peaks, z_trace_later[peaks], "x")
        if PkPos2.size: # if there is a second spike
            PkPrm_Kamino[j,k]=PkProperties2["prominences"][0]
            PkPrm_Kamino_norm[j,k]=PkProperties2["prominences"][0]/PkProperties1["prominences"][0]
        else:
            PkPrm_Kamino[j,k]=0 # if there is no second spike
            PkPrm_Kamino_norm[j,k] = 0
 
           
# plot FC vs. second response amplitude
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Kamino)+1))
label_font_size = 9
trace_width = 2

fig3 = plt.figure(figsize=(6, 6)) 
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)

ax1= fig3.add_subplot(grid[0, 0])

for i in range(len(z0First_space_Kamino)):
    ax1.plot(FC_space_Kamino,PkPrm_Kamino[i,:], linewidth=trace_width,
             color = colors[i], label='Priming cAMP='+str(z0First_space_Kamino[i]))

ax1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
ax1.set_xlabel('Fold change in extracellular cAMP',fontsize=label_font_size)
ax1.set_xscale('log')
ax1.set_title('Second spike prominence')
leg = ax1.legend();

#%% Sgro 2015, w/o noise (sigma=0)
from Sgro2015_agent_and_pop_FUN import Sgro2015_agent

z0First_space_Sgro =np.array([1,2,4,8]) # np.array([0.125,0.25,0.5,1,2,4]) # first period extracellular cAMP level # 
FC_space_Sgro=  np.logspace(0.1, 1.6, num=8) # np.logspace(0.5, 0.5, num=1) 
run_time_space =[0]# np.arange(0,0,1) # run 10 times and plot the mean of each case
# Initialize 

PkPrm_Sgro_woNoise = np.zeros((len(z0First_space_Sgro), len(FC_space_Sgro), len(run_time_space))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus
PkPrm_Sgro_woNoise_norm = np.zeros((len(z0First_space_Sgro), len(FC_space_Sgro), len(run_time_space)))

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'flux_thrs':0}
Na=3.5;  # normalization factor of A
A0=-1.5; R0=-0.5
Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)

dt=0.001 ; t_tot=30*Nt_Sgro; t=list(np.arange(0,t_tot,dt))

# j_test=[2]; k_test=[3]

for j in range(len(z0First_space_Sgro)):
    signal_trace=z0First_space_Sgro[j]*np.zeros(len(t))
    for k in range(len(FC_space_Sgro)):       
        stim_time_step1 = int(round(1/3*t_tot/dt))
        stim_time_step2=int(round(2/3*t_tot/dt)) # at this time second step input is applied
        
        signal_trace[stim_time_step1:stim_time_step2] = z0First_space_Sgro[j]
        signal_trace[stim_time_step2:] = FC_space_Sgro[k]*z0First_space_Sgro[j]
        
        for test in run_time_space:
            A_trace_orig=[A0]; R_trace_orig=[R0]
            for i in range(len(t)-1):
                A_now=A_trace_orig[i]
                R_now=R_trace_orig[i]
                signal_now=signal_trace[i]
        
                A_next,R_next,r_now=Sgro_agent.update(dt,signal_now, r=0)
                A_trace_orig.append(A_next)
                R_trace_orig.append(R_next)
       
            A_trace_offset=1.5
            A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
            A_trace_plot=(A_trace_orig - Nh_Sgro_offset)/Nh_Sgro;
            t_plot_Sgro = np.array(t)/Nt_Sgro
            
#            # check traces
#            fig3 = plt.figure(figsize=(6, 6))
#            grid = plt.GridSpec(2, 1, wspace=0.3, hspace=0.2)
#            ax1= fig3.add_subplot(grid[0, 0])
#            ax1.plot(t_plot_Sgro,signal_trace)
#            ax1.set_ylabel('extracellular cAMP')
#            ax1.set_title('cAMP from '+str(z0First_space_Sgro[j])+' to FC '+ str(FC_space_Sgro[k]))
#            ax2= fig3.add_subplot(grid[1, 0])
#            ax2.plot(t_plot_Sgro,A_trace_plot)
#            ax2.set_ylabel('Activator')
#            ax2.set_ylabel('Time')
#            plt.show()
            
#            # check traces
#            plt.figure(figsize=(4,2.5))
#            # plt.plot(t_plot_Kamino,x_trace,t_plot_Kamino,y_trace,t_plot_Kamino,signal_trace)
#            plt.plot(t_plot_Sgro/Nt_Sgro,A_trace_plot)
#            plt.xlabel('Time')
#            plt.ylabel('A')
#            plt.title('cAMP from  '+str(z0First_space_Sgro[j])+' to '+str(FC_space_Sgro[k]*z0First_space_Sgro[j]))
#            # plt.gca().legend(('x','y','z'))
#            plt.show()
        
            A_trace_first=A_trace_plot[stim_time_step1:stim_time_step2] 
            A_trace_second=A_trace_plot[stim_time_step2:] # the second part of trace, second spike
            PkPos1, PkProperties1 = find_peaks(A_trace_first, prominence=(0.1,2))
            PkPos2, PkProperties2 = find_peaks(A_trace_second, prominence=(0.1,2))
            # Check find_peaks
            # plt.plot(z_trace_later)
            # plt.plot(peaks, z_trace_later[peaks], "x")
            if PkPos2.size: # if there is a second spike
                PkPrm_Sgro_woNoise[j,k,test]=PkProperties2["prominences"][0]
                PkPrm_Sgro_woNoise_norm[j,k,test]=PkProperties2["prominences"][0]/PkProperties1["prominences"][0]
            else:
                PkPrm_Sgro_woNoise[j,k,test]=0 # if there is no second spike
                PkPrm_Sgro_woNoise_norm[j,k,test] = 0
 
PkPrm_Sgro_mean=np.mean(PkPrm_Sgro_woNoise,axis=2)
# plot FC vs. second response amplitude
#fig3 = plt.figure(figsize=(16, 14))
#grid = plt.GridSpec(2, 2, wspace=0.1, hspace=0.35)
label_font_size = 9
trace_width = 2

fig3 = plt.figure(figsize=(6, 4.5))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)

colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Sgro)+1))
title_font_size = 22
label_font_size = 20
tick_font_size = 16
legend_font_size = 15
trace_width = 5
ax3= fig3.add_subplot(grid[0, 0])
for i in range(len(z0First_space_Sgro)):
    ax3.plot(FC_space_Sgro,PkPrm_Sgro_woNoise[i,:], color = colors[i],
             linewidth=trace_width,label='Prime Conc.'+str(z0First_space_Sgro[i]))
ax3.set_xscale('log')
ax3.set_ylim([-0.3,1])
ax3.set_xlim([FC_space_Sgro[0],FC_space_Sgro[-1]])
ax3.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
ax3.set_xlabel(r'$cAMP_{e}$'+' fold change',fontsize=label_font_size)
# ax3.set_xscale('log')
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Sgro 2015, w/o noise', fontdict={'fontsize': title_font_size})
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
leg = ax3.legend();
ax3.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})

fig3.savefig('figS12_FCD_Sgro_wonoise_tight_200604.png', bbox_inches='tight') 

#%% Sgro 2015, with noise
from Sgro2015_agent_and_pop_FUN import Sgro2015_agent

z0First_space_Sgro_noise =np.array([1,2,4,8])#([0.5,1,2]) np.array([0.125,0.25,0.5,1,2,4]) # first period extracellular cAMP level # 
FC_space_Sgro_noise=  np.logspace(0.1, 1.6, num=8) # np.logspace(0.5, 0.5, num=1) 
num_of_runs = 20
run_time_space =np.arange(0,num_of_runs,1) # run 10 times and plot the mean of each case
# Initialize 

PkPrm_Sgro_noise = np.zeros((len(z0First_space_Sgro_noise), len(FC_space_Sgro_noise), len(run_time_space))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus
PkPrm_Sgro_noise_norm = np.zeros((len(z0First_space_Sgro_noise), len(FC_space_Sgro_noise), len(run_time_space))) 

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
Na=3.5;  # normalization factor of A
A0=-1.5; R0= -0.5
Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)

dt=0.005 ; t_tot=30*Nt_Sgro; t=list(np.arange(0,t_tot,dt))

# j_test=[2]; k_test=[3]                                 

for j in range(len(z0First_space_Sgro_noise)):
    signal_trace=np.zeros(len(t))
    for k in range(len(FC_space_Sgro_noise)):
        
        stim_time_step1 = int(round(1/3*t_tot/dt))
        stim_time_step2=int(round(2/3*t_tot/dt)) # at this time second step input is applied
        
        signal_trace[stim_time_step1:stim_time_step2] = z0First_space_Sgro_noise[j]
        signal_trace[stim_time_step2:] = FC_space_Sgro_noise[k]*z0First_space_Sgro_noise[j]
        
        for test in run_time_space:
            A_trace_orig=[A0]; R_trace_orig=[R0]
            for i in range(len(t)-1):
                A_now=A_trace_orig[i]
                R_now=R_trace_orig[i]
                signal_now=signal_trace[i]
        
                A_next,R_next,r_now=Sgro_agent.update(dt,signal_now)
                A_trace_orig.append(A_next)
                R_trace_orig.append(R_next)
       
            A_trace_offset=1.5
            A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
            A_trace_plot=(A_trace_orig-A_trace_offset)/Nh_Sgro;
            t_plot_Sgro_noise = np.array(t)/Nt_Sgro
            
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
    
#            A_trace_first=A_trace_plot[stim_time_step1:stim_time_step2] 
#            A_trace_second=A_trace_plot[stim_time_step2:] # the second part of trace, second spike
            A_trace_first=A_trace_plot[stim_time_step1:int(stim_time_step1+1.5*Nt_Sgro/dt)] 
            A_trace_second=A_trace_plot[stim_time_step2:int(stim_time_step2+1.5*Nt_Sgro/dt)]         
            PkPos1, PkProperties1 = find_peaks(A_trace_first, prominence=(0.4,2))
            PkPos2, PkProperties2 = find_peaks(A_trace_second, prominence=(0.4,2))
            # Check find_peaks
            # plt.plot(z_trace_later)
            # plt.plot(peaks, z_trace_later[peaks], "x")
            if PkPos2.size: # if there is a second spike
                PkPrm_Sgro_noise[j,k,test]=PkProperties2["prominences"][0]
                PkPrm_Sgro_noise_norm[j,k,test]=PkProperties2["prominences"][0]/PkProperties1["prominences"][0]
            else:
                PkPrm_Sgro_noise[j,k,test]=0 # if there is no second spike
                PkPrm_Sgro_noise_norm[j,k,test]=0 
    
    print('j='+str(j)+'is done')
    
PkPrm_Sgro_mean_noise=np.mean(PkPrm_Sgro_noise,axis=2)
PkPrm_Sgro_se_noise = np.std(PkPrm_Sgro_noise,axis=2)/sqrt(num_of_runs)
PkPrm_Sgro_std_noise = np.std(PkPrm_Sgro_noise,axis=2)

PkPrm_Sgro_mean_noise_norm=np.mean(PkPrm_Sgro_noise_norm,axis=2)
PkPrm_Sgro_se_noise_norm = np.std(PkPrm_Sgro_noise_norm,axis=2)/sqrt(num_of_runs)
PkPrm_Sgro_std_noise_norm = np.std(PkPrm_Sgro_noise_norm,axis=2)

# plot FC vs. second response amplitude
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Sgro_noise)+1))
fig3 = plt.figure(figsize=(8, 7))
grid = plt.GridSpec(1, 1, wspace=0.1, hspace=0.35)

title_font_size = 20 
label_font_size = 20
tick_font_size = 16
legend_font_size = 12
trace_width = 5
markers=['+', 'd', '2', 'x']
# plot FC vs. second response amplitude, scatter plot
ax3= fig3.add_subplot(grid[0, 0])
for i in range(len(z0First_space_Sgro_noise)):
    ax3.plot(FC_space_Sgro_noise, PkPrm_Sgro_mean_noise[i,:],'o-', color = colors[i], lw = 3, ms = 10,  label='Prime cAMP='+str(z0First_space_Sgro_noise[i]))
    ax3.errorbar(FC_space_Sgro_noise, PkPrm_Sgro_mean_noise[i,:], yerr=PkPrm_Sgro_se_noise[i,:],
                 fmt = 'o', color=colors[i], ecolor= colors[i], elinewidth=3, capsize=10, capthick=3)
for i in range(len(z0First_space_Sgro_noise)):
    z0_now = z0First_space_Sgro_noise[i]
    for j in range(len(FC_space_Sgro_noise)):
        FC_now = FC_space_Sgro_noise[j]
        x = FC_now*np.ones(num_of_runs)
        ax3.plot(x, PkPrm_Sgro_noise[i,j,:],markers[i], color = colors[i], ms = 5)
    
ax3.set_ylim([-0.1,1.3])
ax3.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
ax3.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax3.set_xscale('log')
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Sgro 2015 (with noise)', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
#leg = ax3.legend();
ax3.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})

#%% violin plot
colors = plt.cm.Greens(np.linspace(0,1,len(z0First_space_Sgro_noise)+2))

title_font_size = 20
label_font_size = 18
tick_font_size = 14
legend_font_size = 16
trace_width = 5
markers=['+', 'd', '2', 'x']
# plot FC vs. second response amplitude, violoin plot
fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(cm2inch(50),cm2inch(12)))
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
    ax[i].set_xticklabels(np.around(FC_space_Sgro_noise,decimals=0),fontsize=tick_font_size)
    ax[i].tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax[i].set_title('Prime concentration '+str(round(z0_now,2)), fontdict={'fontsize': label_font_size, 'fontweight': 'medium'})
fig.text(0.5, 0.97, 'Sgro 2015 (with noise)', ha='center', va='center',fontsize=title_font_size,color=sgrocolor,)
fig.text(0.5, 0.02, r'$cAMP_{e}$'+' fold change', ha='center', va='center',fontsize=label_font_size)
fig.text(0.08, 0.5, 'Second peak prominence, A.U.', ha='center', va='center', rotation='vertical',fontsize=label_font_size)


#plt.setp(violin_parts['bodies'], facecolor=colors[i],edgecolor=colors[i],alpha=0.2)
#plt.setp(violin_parts['bodies'], facecolor=colors[i],edgecolor=colors[i])



        # sns.violinplot(ax=ax3, x = x_now, y=PkPrm_Sgro_noise[i,j,:])
        
#%%       
data = [sorted(np.random.normal(0, std, 100)) for std in range(1, 5)]
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(9, 4), sharey=True)

ax1.set_title('Default violin plot')
ax1.set_ylabel('Observed values')
ax1.violinplot(data)
#%%
ax3.set_ylim([-0.1,1.3])
ax3.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
ax3.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax3.set_xscale('log')
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Sgro 2015 (with noise)', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
#leg = ax3.legend();
ax3.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})

#%% Goldbeter 1987 FCD
from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_agent_3var

z0First_space_Gold =np.array([0.1,0.2,0.4,0.8]) #([0.1,0.2,0.4,0.8]) # first period extracellular cAMP level # 
FC_space_Gold=  np.logspace(0.1, 1.6, num=12) # np.logspace(0.5, 0.5, num=1) 

# Initialize 
PkPrm_Gold = np.zeros((len(z0First_space_Gold), len(FC_space_Gold))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus
PkPrm_Gold_norm = np.zeros((len(z0First_space_Gold), len(FC_space_Gold)))

k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10; L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01; theta=0.01
e= 1 # compared to 0.108
q=4000
sig= 0.6# compared to 0.57
v=12; k= 4 # k prime in the paper
ki=1.7 # compared to 0.958
kt=0.9
kc=5.4 # compared to 3.58
h=5

Goldbeter3AgentParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}

p0=0.8; a0=3; b0=0.9; g0=0

dt=0.0005; t_tot=30*Nt_Goldbeter; t=list(np.arange(0,t_tot,dt))

j_test=[1]; k_test=[1,3,5]

for j in range(len(z0First_space_Gold)): # j_test: # 
    signal_trace = np.zeros(len(t)) 
    for k in range(len(FC_space_Gold)): # k_test: # 
        stim_time_step1 = int(round(1/3*t_tot/dt))
        stim_time_step2=int(round(2/3*t_tot/dt)) # at this time second step input is applied
        
        signal_trace[stim_time_step1:stim_time_step2] = z0First_space_Gold[j]
        signal_trace[stim_time_step2:] = FC_space_Gold[k]*z0First_space_Gold[j]
        
        p_trace=[p0]; b_trace=[b0]; g_trace=[g0]
        Goldbeter3_agent=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
        for i in range(len(t)-1):
            p_now=p_trace[i]
            p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,signal_trace[i])
            p_trace.append(p_next)
            b_trace.append(b_next)
            g_trace.append(g_next)
   
        b_trace = np.array(b_trace)/Nh_Goldbeter; # b_trace = b_trace/np.amax(b_trace)
        t_plot_Gold = np.array(t)/Nt_Goldbeter
        
        
        b_trace_first=b_trace[stim_time_step1:stim_time_step2]
        b_trace_second=b_trace[stim_time_step2:] # the second part of trace, second spike
        PkPos1, PkProperties1 = find_peaks(b_trace_first, prominence=(0,3000/Nh_Goldbeter))
        PkPos2, PkProperties2 = find_peaks(b_trace_second, prominence=(0,3000/Nh_Goldbeter))
        
        Pk1Pos = PkPos1[0]+stim_time_step1; Pk2Pos = PkPos2[0] + stim_time_step2
        Pk1max = b_trace[Pk1Pos]; Pk1min = b_trace[Pk1Pos] -PkProperties1["prominences"][0]
        Pk2max = b_trace[Pk2Pos]; Pk2min = b_trace[Pk2Pos] -PkProperties2["prominences"][0]
#        # check traces & peaks
#        fig3 = plt.figure(figsize=(4, 3))
#        grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
#        ax1= fig3.add_subplot(grid[0, 0])
#        ax1.plot(t_plot_Gold,signal_trace)
#        ax1.set_ylabel('extracellular cAMP')
#        ax1.set_title('cAMP from '+str(z0First_space_Gold[j])+' to FC '+ str(FC_space_Gold[k]))
#        ax2= fig3.add_subplot(grid[1:, 0])
#        ax2.plot(t_plot_Gold,b_trace)
#        ax2.plot([Pk1Pos*dt/Nt_Goldbeter,Pk1Pos*dt/Nt_Goldbeter],[Pk1min,Pk1max],color = 'g',linewidth=2.5)
#        ax2.text(Pk1Pos*dt/Nt_Goldbeter, 5, str(round(Pk1max-Pk1min,2)) , rotation=90, va='center')
#        ax2.plot([Pk2Pos*dt/Nt_Goldbeter,Pk2Pos*dt/Nt_Goldbeter],[Pk2min,Pk2max],color = 'g',linewidth=2.5)
#        ax2.text(Pk2Pos*dt/Nt_Goldbeter, 5, str(round(Pk2max-Pk2min,2)) , rotation=90, va='center')
#        
#        ax2.set_ylabel('beta, [cAMP]cyto')
#        ax2.set_ylabel('Time')
#        plt.show()

        if PkPos2.size: # if there is a second spike
            PkPrm_Gold_norm[j,k]=PkProperties2["prominences"][0]/PkProperties1["prominences"][0]
            PkPrm_Gold[j,k]=PkProperties2["prominences"][0]
        else:
            PkPrm_Gold_norm[j,k]= 0 # if there is no second spike
            PkPrm_Gold[j,k]= 0 
#        if  FC_space_Gold[k] <= z0First_space_Gold[j]:
#            PkPrm_Gold[j,k] = 0            

# PkPrm_Gold_mean=np.mean(PkPrm_Gold,axis=2)
# plot FC vs. second response amplitude
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Gold)+1))
label_font_size = 15
trace_width = 2
tick_font_size = 14

fig4 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)

ax1= fig4.add_subplot(grid[0, 0])

for i in range(len(z0First_space_Gold)):
    ax1.plot(FC_space_Gold,PkPrm_Gold[i,:], linewidth=trace_width, color=colors[i],
             label='Priming cAMP='+str(z0First_space_Gold[i]))

ax1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
ax1.set_xscale('log')
#ax1.set_xlabel('Second step of extracellular cAMP',fontsize=label_font_size)
ax1.set_xlabel('Second step of extracellular cAMP',fontsize=label_font_size)
ax1.set_title('Goldbeter 1987') # , 2nd peak width normalized
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
leg = ax1.legend();      

#%% Goldbeter 1987 second stim concentration
from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_agent_3var

z0First_space_Gold =np.array([0.1,0.2,0.4,0.8]) #([0.1,0.2,0.4,0.8]) # first period extracellular cAMP level # 
scd_space_Gold=  np.logspace(-0.5, 1.6, num=8) 

# Initialize 
PkPrm_Gold = np.zeros((len(z0First_space_Gold), len(scd_space_Gold))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus
PkPrm_Gold_norm = np.zeros((len(z0First_space_Gold), len(scd_space_Gold)))

k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10; L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01; theta=0.01
e= 1 # compared to 0.108
q=4000
sig= 0.6# compared to 0.57
v=12; k= 4 # k prime in the paper
ki=1.7 # compared to 0.958
kt=0.9
kc=5.4 # compared to 3.58
h=5

Goldbeter3AgentParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}

p0=0.8; a0=3; b0=0.9; g0=0

dt=0.0005; t_tot=30*Nt_Goldbeter; t=list(np.arange(0,t_tot,dt))

# j_test=[1]; k_test=[1,3,5]

# show traces of these parameter indexes, 1st column z0_First, 2nd column FC
show_trace_idx = np.array([[3,1],[3,3],[3,5],[3,7]])

for j in range(len(z0First_space_Gold)): # j_test: # 
    signal_trace = np.zeros(len(t)) 
    for k in range(len(scd_space_Gold)): # k_test: # 
        stim_time_step1 = int(round(1/3*t_tot/dt))
        stim_time_step2=int(round(2/3*t_tot/dt)) # at this time second step input is applied
        
        signal_trace[stim_time_step1:stim_time_step2] = z0First_space_Gold[j]
        signal_trace[stim_time_step2:] = scd_space_Gold[k]
        
        p_trace=[p0]; b_trace=[b0]; g_trace=[g0]
        Goldbeter3_agent=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
        for i in range(len(t)-1):
            p_now=p_trace[i]
            p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,signal_trace[i])
            p_trace.append(p_next)
            b_trace.append(b_next)
            g_trace.append(g_next)
   
        b_trace = np.array(b_trace)/Nh_Goldbeter; # b_trace = b_trace/np.amax(b_trace)
        t_plot_Gold = np.array(t)/Nt_Goldbeter
        
        
        b_trace_first=b_trace[stim_time_step1:stim_time_step2]
        b_trace_second=b_trace[stim_time_step2:] # the second part of trace, second spike
        PkPos1, PkProperties1 = find_peaks(b_trace_first, prominence=(0,300))
        PkPos2, PkProperties2 = find_peaks(b_trace_second, prominence=(0,300))
        
        Pk1Pos = PkPos1[0]+stim_time_step1; 
        Pk1max = b_trace[Pk1Pos]; Pk1min = b_trace[Pk1Pos] -PkProperties1["prominences"][0]
        if PkPos2.size!=0:
            Pk2Pos = PkPos2[0] + stim_time_step2
            Pk2max = b_trace[Pk2Pos]; Pk2min = b_trace[Pk2Pos] -PkProperties2["prominences"][0]
        
        # show and check selected traces & peaks
        if j in show_trace_idx[:,0] and k in show_trace_idx[:,1]:       
            fig3 = plt.figure(figsize=(6, 4))
            grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
            ax1= fig3.add_subplot(grid[0, 0])
            ax1.plot(t_plot_Gold,signal_trace)
            ax1.set_ylabel('extracellular cAMP')
            ax1.set_title('cAMP from '+str(z0First_space_Gold[j])+' to '+ str(scd_space_Gold[k]))
            ax2= fig3.add_subplot(grid[1:, 0])
            ax2.plot(t_plot_Gold,b_trace)
            ax2.plot([Pk1Pos*dt/Nt_Goldbeter,Pk1Pos*dt/Nt_Goldbeter],[Pk1min,Pk1max],color = 'g',linewidth=2.5)
            ax2.text(Pk1Pos*dt/Nt_Goldbeter, 5, str(round(Pk1max-Pk1min,2)) , rotation=90, va='center')
            if size(PkPos2)!=0:
                ax2.plot([Pk2Pos*dt/Nt_Goldbeter,Pk2Pos*dt/Nt_Goldbeter],[Pk2min,Pk2max],color = 'g',linewidth=2.5)
                ax2.text(Pk2Pos*dt/Nt_Goldbeter, 5, str(round(Pk2max-Pk2min,2)) , rotation=90, va='center')
            ax2.set_ylabel(r'[cAMP]i')
            ax2.set_ylabel('Time, A.U.')
            plt.show()

        if PkPos2.size: # if there is a second spike
            PkPrm_Gold_norm[j,k]=PkProperties2["prominences"][0]/PkProperties1["prominences"][0]
            PkPrm_Gold[j,k]=PkProperties2["prominences"][0]
        else:
            PkPrm_Gold_norm[j,k]= 0 # if there is no second spike
            PkPrm_Gold[j,k]= 0 
#        if  FC_space_Gold[k] <= z0First_space_Gold[j]:
#            PkPrm_Gold[j,k] = 0            

# PkPrm_Gold_mean=np.mean(PkPrm_Gold,axis=2)

# plot scd stimulation concentration vs. second response amplitude
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Gold)+1))
trace_width = 2
title_font_size = 18
label_font_size = 16
tick_font_size = 14

fig3 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)

ax1= fig3.add_subplot(grid[0, 0])

for i in range(len(z0First_space_Gold)):
    
    ax1.plot(scd_space_Gold,PkPrm_Gold[i,:], linewidth=trace_width, color = colors[i],
             label='Priming cAMP='+str(z0First_space_Gold[i]))

ax1.set_ylabel( 'second spike prominence',fontsize=label_font_size)
ax1.set_xscale('log')
#ax1.set_xlabel('Fold change in extracellular cAMP',fontsize=label_font_size)
ax1.set_xlim([scd_space_Gold[0],scd_space_Gold[-1]])
ax1.set_xlabel(r'Second step [$cAMP_{e}$]',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Martiel 1987', color = 'k',fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax1.legend();   

#%% Maeda and Loomis 2004 FCD
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


z0First_space_Loomis =np.array([0.1,0.2,0.4,0.8]) #np.array([0.5,10]) # first period extracellular cAMP level # 
FC_space_Loomis=  np.logspace(0.1, 1.6, num=8) # np.logspace(0.5, 0.5, num=1) 
#z0Second_space_Loomis = np.logspace(-0.7, 1.7, num=8)

# Initialize 
PkPrm_Loomis = np.zeros((len(z0First_space_Loomis), len(FC_space_Loomis))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus
PkPrm_Loomis_norm = np.zeros((len(z0First_space_Loomis), len(FC_space_Loomis)))

dt=0.0005; t_tot=30*Nt_Maeda; t=list(np.arange(0,t_tot,dt))
show_trace_idx = np.array([[3,1],[3,3],[3,5],[3,7]])
#  k_test=np.array([7])
# j_test=[1];k_test=[0,1,2,3,4,5,6,7]

for j in range(len(z0First_space_Loomis)):  # j_test:#
    signal_trace=z0First_space_Loomis[j]*np.ones(len(t))
    for k in range(len(FC_space_Loomis)):# k_test: #
        stim_time_step1 = int(round(1/3*t_tot/dt))
        stim_time_step2=int(round(2/3*t_tot/dt)) # at this time second step input is applied
        
        signal_trace[stim_time_step1:stim_time_step2] = z0First_space_Loomis[j]
        signal_trace[stim_time_step2:] = FC_space_Loomis[k]*z0First_space_Loomis[j]
        
        ACA_trace=[ACA0]; PKA_trace=[PKA0]; ERK2_trace=[ERK20]; RegA_trace=[RegA0]; 
        cAMPi_trace=[cAMPi0]; cAMPe_trace=[cAMPe0]; CAR1_trace=[CAR10]
        
        Maeda_agent=MaedaLoomis2004_agent([1,1],state0,MaedaAgentParam)
        
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
        t_plot_Loomis = np.array(t)/Nt_Maeda
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

        
        cAMPi_trace_first=cAMPi_trace[stim_time_step1:stim_time_step2]
        cAMPi_trace_second=cAMPi_trace[stim_time_step2:] # the second part of trace, second spike
        PkPos1, PkProperties1 = find_peaks(cAMPi_trace_first, prominence=(0,100))
        PkPos2, PkProperties2 = find_peaks(cAMPi_trace_second, prominence=(0,100))
        
        Pk1Pos = PkPos1[0]+stim_time_step1; Pk2Pos = PkPos2[0] + stim_time_step2
        Pk1max = cAMPi_trace[Pk1Pos]; Pk1min = cAMPi_trace[Pk1Pos] -PkProperties1["prominences"][0]
        Pk2max = cAMPi_trace[Pk2Pos]; Pk2min = cAMPi_trace[Pk2Pos] -PkProperties2["prominences"][0]
        
        if j in show_trace_idx[:,0] and k in show_trace_idx[:,1]:              
            # Check find_peaks
            fig3 = plt.figure(figsize=(4, 4))
            grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
            ax= fig3.add_subplot(grid[0, 0])
            ax.plot(cAMPi_trace_second)
            ax.plot(PkPos2, cAMPi_trace_second[PkPos2], "x")
            # ax.set_ylim([0.5,25])
            ax.set_title('cAMP from '+str(z0First_space_Loomis[j])+' to FC '+ str(FC_space_Loomis[k]))
            # ax.set_title('cAMP from '+str(z0First_space_Loomis[j])+' to '+ str(z0Second_space_Loomis[k]))
            plt.show()
        
        if PkPos2.size: # if there is a second spike
            PkPrm_Loomis[j,k]=PkProperties2["prominences"][0]
            PkPrm_Loomis_norm[j,k]=PkProperties2["prominences"][0]/PkProperties1["prominences"][0]
        else:
            PkPrm_Loomis[j,k]=0 # if there is no second spike
            PkPrm_Loomis_norm[j,k]=0
        
                

# Plot the results
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Gold)+1))
trace_width = 2
title_font_size = 18
label_font_size = 16
tick_font_size = 14

fig3 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)

ax1= fig3.add_subplot(grid[0, 0])

for i in range(len(z0First_space_Loomis)):
    
    ax1.plot(FC_space_Loomis,PkPrm_Loomis[i,:], linewidth=trace_width, color = colors[i],
             label='Priming cAMP='+str(z0First_space_Loomis[i]))

ax1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
ax1.set_xscale('log')
ax1.set_xlim([FC_space_Loomis[0],FC_space_Loomis[-1]])
ax1.set_xlabel(r'Fold change in [$cAMP_{e}$]',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Maeda 2004', color = 'k',fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax1.legend();


#%% Maeda and Loomis 2004 second stim concentration
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


z0First_space_Loomis =np.array([0.1,0.2,0.4,0.8]) #np.array([0.5,10]) # first period extracellular cAMP level # 
scd_space_Loomis=  np.logspace(-0.5, 1.6, num=8) # np.logspace(0.5, 0.5, num=1) 


# Initialize 
PkPrm_Loomis = np.zeros((len(z0First_space_Loomis), len(scd_space_Loomis))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus
PkPrm_Loomis_norm = np.zeros((len(z0First_space_Loomis), len(scd_space_Loomis)))

dt=0.0005; t_tot=30*Nt_Maeda; t=list(np.arange(0,t_tot,dt))

#  k_test=np.array([7])
# j_test=[1];k_test=[0,1,2,3,4,5,6,7]
# show traces of these parameter indexes, 1st column z0_First, 2nd column FC
show_trace_idx = np.array([[3,1],[3,3],[3,5],[3,7]])

for j in range(len(z0First_space_Loomis)):  # j_test:#
    signal_trace=z0First_space_Loomis[j]*np.ones(len(t))
    for k in range(len(scd_space_Loomis)):# k_test: #
        stim_time_step1 = int(round(1/3*t_tot/dt))
        stim_time_step2=int(round(2/3*t_tot/dt)) # at this time second step input is applied
        
        signal_trace[stim_time_step1:stim_time_step2] = z0First_space_Loomis[j]
        signal_trace[stim_time_step2:] = scd_space_Loomis[k]
        
        ACA_trace=[ACA0]; PKA_trace=[PKA0]; ERK2_trace=[ERK20]; RegA_trace=[RegA0]; 
        cAMPi_trace=[cAMPi0]; cAMPe_trace=[cAMPe0]; CAR1_trace=[CAR10]
        
        Maeda_agent=MaedaLoomis2004_agent([1,1],state0,MaedaAgentParam)
        
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
        t_plot_Loomis = np.array(t)/Nt_Maeda
        
#            # check traces
#            fig3 = plt.figure(figsize=(3, 3))
#            grid = plt.GridSpec(2, 1, wspace=0.3, hspace=0.2)
#            ax1= fig3.add_subplot(grid[0, 0])
#            ax1.plot(t_plot_Loomis,signal_trace)
#            ax1.set_ylabel('extracellular cAMP')
#            # ax1.set_title('cAMP from '+str(z0First_space_Loomis[j])+' to '+ str(scd_space_Loomis[k]))
#            ax1.set_title('cAMP from '+str(z0First_space_Loomis[j])+' to secondary concentration '+ str(scd_space_Loomis[k]))
#            ax2= fig3.add_subplot(grid[1, 0])
#            ax2.plot(t_plot_Loomis,cAMPi_trace)
#            ax2.set_ylim([0,25])
#            ax2.set_ylabel(' [cAMP]i')
#            ax2.set_ylabel('Time, A.U.')
#            plt.show()

        
        cAMPi_trace_first=cAMPi_trace[stim_time_step1:stim_time_step2]
        cAMPi_trace_second=cAMPi_trace[stim_time_step2:] # the second part of trace, second spike
        PkPos1, PkProperties1 = find_peaks(cAMPi_trace_first, prominence=(0,100))
        PkPos2, PkProperties2 = find_peaks(cAMPi_trace_second, prominence=(0,100))
        
        Pk1Pos = PkPos1[0]+stim_time_step1; Pk2Pos = PkPos2[0] + stim_time_step2
        Pk1max = cAMPi_trace[Pk1Pos]; Pk1min = cAMPi_trace[Pk1Pos] -PkProperties1["prominences"][0]
        Pk2max = cAMPi_trace[Pk2Pos]; Pk2min = cAMPi_trace[Pk2Pos] -PkProperties2["prominences"][0]
        
        if j in show_trace_idx[:,0] and k in show_trace_idx[:,1]:       
            # Check find_peaks
            fig3 = plt.figure(figsize=(4, 4))
            grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
            ax= fig3.add_subplot(grid[0, 0])
            ax.plot(cAMPi_trace_second)
            ax.plot(PkPos2, cAMPi_trace_second[PkPos2], "x")
            ax.set_ylim([0.5,25])
            ax.set_title('cAMP from '+str(z0First_space_Loomis[j])+' to '+ str(scd_space_Loomis[k]))
            # ax.set_title('cAMP from '+str(z0First_space_Loomis[j])+' to '+ str(z0Second_space_Loomis[k]))
            plt.show()
        
        if PkPos2.size: # if there is a second spike
            PkPrm_Loomis[j,k]=PkProperties2["prominences"][0]
            PkPrm_Loomis_norm[j,k]=PkProperties2["prominences"][0]/PkProperties1["prominences"][0]
        else:
            PkPrm_Loomis[j,k]=0 # if there is no second spike
            PkPrm_Loomis_norm[j,k]=0
        

# plot scd stimulation concentration vs. second response amplitude
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Gold)+1))
trace_width = 2
title_font_size = 18
label_font_size = 16
tick_font_size = 14

fig3 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)

ax1= fig3.add_subplot(grid[0, 0])

for i in range(len(z0First_space_Loomis)):
    
    ax1.plot(scd_space_Loomis,PkPrm_Loomis[i,:], linewidth=trace_width, color = colors[i],
             label='Priming cAMP='+str(z0First_space_Loomis[i]))

ax1.set_ylabel( 'second spike prominence',fontsize=label_font_size)
ax1.set_xscale('log')
#ax1.set_xlabel('Fold change in extracellular cAMP',fontsize=label_font_size)
ax1.set_xlim([scd_space_Loomis[0],scd_space_Loomis[-1]])
ax1.set_xlabel(r'Second step [$cAMP_{e}$]',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Maeda 2004', color = 'k',fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax1.legend();

#%% save outputs as npz files
np.savez('figS5_sc_FCD_2ndconc_042320.npz', 
         z0First_space_Gold = z0First_space_Gold , 
         scd_space_Gold= scd_space_Gold,
         PkPrm_Gold = PkPrm_Gold,
         z0First_space_Loomis = z0First_space_Loomis,
         scd_space_Loomis = scd_space_Loomis,
         PkPrm_Loomis = PkPrm_Loomis)
#%% Plot all 4 models
#title_font_size = 20 
#label_font_size = 20
#tick_font_size = 16
#legend_font_size = 15
#trace_width = 5
#
##fig3 = plt.figure(figsize=(16, 14))
##grid = plt.GridSpec(2, 2, wspace=0.15, hspace=0.25)
#fig3 = plt.figure(figsize=(12,13))
#grid = plt.GridSpec(2, 2, wspace=0.1, hspace=0.3)
#
#ax1= fig3.add_subplot(grid[0, 0])
#for i in range(len(z0First_space_Gold)):
#    ax1.plot(FC_space_Gold,PkPrm_Gold[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Gold[i]))
#
#ax1.set_ylim([-0.3,1])
##ax1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
##ax1.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
#ax1.set_xscale('log')
#ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
#ax1.set_title('Goldbeter 1987', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
#leg = ax1.legend();
#ax1.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})
#
#ax2= fig3.add_subplot(grid[0, 1])
#for i in range(len(z0First_space_Loomis)):
#    ax2.plot(FC_space_Loomis,PkPrm_Loomis[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Loomis[i]))
#
#ax2.set_ylim([-1,3])
##ax2.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
##ax2.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
#ax2.set_xscale('log')
#ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
#ax2.set_title('Laub & Loomis 1998', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
#leg = ax2.legend();
#ax2.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})
#
#
#ax3= fig3.add_subplot(grid[1, 0])
#for i in range(len(z0First_space_Sgro)):
#    ax3.plot(FC_space_Sgro_noise,PkPrm_Sgro_mean_noise[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Sgro_noise[i]))
#
#ax3.set_ylim([-0.3,0.9])
##ax3.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
##ax3.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
#ax3.set_xscale('log')
#ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
#ax3.set_title('Sgro 2015 (with noise)', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
#leg = ax3.legend();
#ax3.legend( frameon=False,loc='bottom',ncol=2,prop={'size': legend_font_size})
#
#
#ax4= fig3.add_subplot(grid[1, 1])
#for i in range(len(z0First_space_Kamino)):
#    ax4.plot(FC_space_Kamino,PkPrm_Kamino[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Kamino[i]))
#ax4.set_ylim([-0.08,0.27])    
##ax4.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
##ax4.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
#ax4.set_xscale('log')
#ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
#ax4.set_title('Kamino & Sawai 2017', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
#leg = ax4.legend();
#ax4.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})
#
#fig3.text(0.5, 0.04, r'$cAMP_{ext}$'+' fold change', ha='center', va='center',fontsize=label_font_size)
#fig3.text(0.06, 0.5, 'Second spike prominence', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
#%% Load experimental data
my_dir = r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/'
Kamino_FCD = pd.read_excel(my_dir+r'Kamino_FCD_exp_data.xlsx',sheet_name='Sheet1')

#%% load saved npz output file
npzfile = np.load('single_cell_FCD_200322.npz')
FC_space_Gold = npzfile['FC_space_Gold']
z0First_space_Gold = npzfile['z0First_space_Gold']
PkPrm_Gold = npzfile['PkPrm_Gold']
PkPrm_Gold_norm = npzfile['PkPrm_Gold_norm']

FC_space_Loomis = npzfile['FC_space_Loomis']
z0First_space_Loomis = npzfile['z0First_space_Loomis']
PkPrm_Loomis= npzfile['PkPrm_Loomis']
PkPrm_Loomis_norm= npzfile['PkPrm_Loomis_norm']

FC_space_Sgro = npzfile['FC_space_Sgro']
z0First_space_Sgro = npzfile['z0First_space_Sgro']
PkPrm_Sgro_woNoise = npzfile['PkPrm_Sgro_woNoise']
PkPrm_Sgro_woNoise_norm = npzfile['PkPrm_Sgro_woNoise_norm']
FC_space_Sgro_noise = npzfile['FC_space_Sgro_noise']
z0First_space_Sgro_noise = npzfile['z0First_space_Sgro_noise']
PkPrm_Sgro_noise = npzfile['PkPrm_Sgro_noise']
PkPrm_Sgro_noise_norm = npzfile['PkPrm_Sgro_noise_norm']
PkPrm_Sgro_mean_noise = npzfile['PkPrm_Sgro_mean_noise']
PkPrm_Sgro_se_noise = npzfile['PkPrm_Sgro_se_noise']
FC_space_Kamino = npzfile['FC_space_Kamino']
z0First_space_Kamino = npzfile['z0First_space_Kamino']
PkPrm_Kamino = npzfile['PkPrm_Kamino']
PkPrm_Kamino_norm = npzfile['PkPrm_Kamino_norm']

#%% Panel B: demonstration of concepts
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

# check traces
plt.figure()
# plt.plot(t_plot_Kamino,x_trace,t_plot_Kamino,y_trace,t_plot_Kamino,signal_trace)
plt.plot(t_plot_Kamino_B,y_trace_B)
plt.xlabel('Time')
plt.ylabel('x,y,z')
plt.title('cAMP from  '+str(z0First)+' to '+str(FC*z0First))
# plt.gca().legend(('x','y','z'))
plt.show()

#%% Plot all 4 models
abcd_font_size = 28
title_font_size = 20
label_font_size = 18 # 16
sublabel_font_size = 18
tick_font_size = 16
legend_font_size = 14
trace_width = 3

#abcd_font_size = 28
#label_font_size=24
#title_font_size = 26
#sublabel_font_size = 22
#trace_width=3
#tick_font_size=20

mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00'] 
          
#fig3 = plt.figure(figsize=(16, 14))
#grid = plt.GridSpec(2, 2, wspace=0.15, hspace=0.25)
fig3 = plt.figure(figsize=(9,18))
grid = plt.GridSpec(4, 2, wspace=0.3, hspace=0.9)

#ax0= fig3.add_subplot(grid[0, 0])
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
ax0.set_xlabel('Fold Change in [cAMP]',fontsize=label_font_size)
ax0.set_xscale('log')
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_title('Experiment', color = 'k',fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
leg = ax0.legend();
ax0.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})
ax0.text(-0.1, 1.15, 'A', ha='center',va='center',
         transform = ax0.transAxes, color = 'b', fontsize=abcd_font_size)
ax0.text(-0.27, 0.5, 'Response \n Amplitude, A.U.', ha='center',va='center', rotation ='vertical',
         transform = ax0.transAxes, fontsize=label_font_size)

ax0u = fig3.add_axes([0.12, 0.88, 0.305,0.075],xticks=[0,10,20,30],yticks=[0,1,2,3])
ax0u.plot(t_plot_Kamino_B,signal_trace, 'k',linewidth=trace_width)
ax0u.set_xlim([0,30]); ax0u.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0u.text(-0.25 , 1.3, 'B', ha='center',va='center',
     transform = ax0u.transAxes, color = 'b', fontsize=abcd_font_size)
ax0u.text(-0.28, 0.5, r'$cAMP_{e}$'+'\n input',ha='center',va='center', rotation='vertical',
     transform = ax0u.transAxes, color = 'k', fontsize=sublabel_font_size)

ax0l = fig3.add_axes([0.12, 0.77, 0.305,0.075],xticks=[0,10,20,30],yticks=[0,0.1,0.2])
ax0l.plot(t_plot_Kamino_B, y_trace_B, linewidth=trace_width)
ax0l.set_xlim([0,30]); ax0l.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0l.set_xlabel('Time, A.U.', size=sublabel_font_size)
ax0l.text(-0.28, 0.5, r'$cAMP_{i}$'+'\n response',ha='center',va='center', rotation='vertical',
     transform = ax0l.transAxes, color = 'k', fontsize=sublabel_font_size)
#ax0l.annotate('annotate', xy=(0.1,0.01), xytext=(1,1),rotation = 'vertical',
#            arrowprops={'arrowstyle': '<->'}, va='center')

ax1= fig3.add_subplot(grid[1, 0])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Gold)+1) )
for i in range(len(z0First_space_Gold)):
    ax1.plot(FC_space_Gold,PkPrm_Gold[i,:], color=colors[i],
             linewidth=trace_width,label='Priming Conc. '+str(int(z0First_space_Gold[i]*10))+' , A.U.')
ax1.set_ylim([0,250]);# ax1.set_xlim([0,10**1.5])
#ax1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#ax1.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax1.set_xscale('log')
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Martiel 1987', color=mycolors[0],
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax1.text(-0.1, 1.15, 'C', ha='center',va='center',
         transform = ax1.transAxes, color = 'g', fontsize=abcd_font_size)
#leg = ax1.legend();
#ax1.legend( frameon=False,loc='bottom center',ncol=1,prop={'size': legend_font_size})

ax2= fig3.add_subplot(grid[1, 1])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Loomis)+1) )
for i in range(len(z0First_space_Loomis)):
    ax2.plot(FC_space_Loomis,PkPrm_Loomis[i,:],color=colors[i], linewidth=trace_width,label='Prime Conc.'+str(z0First_space_Loomis[i]))

ax2.set_ylim([0,2]); # ax2.set_xlim([0,10**1.5])
#ax2.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#ax2.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax2.set_xscale('log')
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Maeda 2004', color=mycolors[1],
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax2.text(-0.1, 1.18, 'D', ha='center',va='center',
         transform = ax2.transAxes, color = 'g', fontsize=abcd_font_size)
#leg = ax2.legend();
#ax2.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})

# Sgro w/o noise
ax3= fig3.add_subplot(grid[3, 0])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Sgro)+1) )
for i in range(len(z0First_space_Sgro)):
    ax3.plot(FC_space_Sgro_noise,PkPrm_Sgro_woNoise[i,:],color=colors[i], linewidth=trace_width,label='Prime Conc.'+str(z0First_space_Sgro_noise[i]))

ax3.set_ylim([-0.2,1]); # ax3.set_xlim([0,10**1.5])
#ax3.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
#ax3.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
ax3.set_xscale('log')
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Sgro 2015 (w/o noise)',  color=mycolors[5],
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax3.text(-0.1, 1.18, 'G', ha='center',va='center',
         transform = ax3.transAxes, color = 'g', fontsize=abcd_font_size)
#leg = ax3.legend();
#ax3.legend( frameon=False,loc='bottom',ncol=2,prop={'size': legend_font_size})

# Sgro with noise
ax5= fig3.add_subplot(grid[2, 0])
for i in range(len(z0First_space_Sgro_noise)):
    ax5.plot(FC_space_Sgro_noise, PkPrm_Sgro_noise_norm[i,:],'o-', color = colors[i], lw = trace_width, 
             ms = 4,  label='Prime Conc.'+str(z0First_space_Sgro_noise[i]))
    ax5.errorbar(FC_space_Sgro_noise, PkPrm_Sgro_noise_norm[i,:], yerr=PkPrm_Sgro_se_noise[i,:],
                 fmt = 'o', color=colors[i], ecolor= colors[i], elinewidth=trace_width, capsize=5, capthick=2)
ax5.set_ylim([-0.2,1]);# ax5.set_xlim([0,10**1.5])
ax5.set_xscale('log')
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5.set_title('Sgro 2015 (noise)', color=mycolors[5],
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax5.text(-0.1, 1.18, 'E', ha='center',va='center',
         transform = ax5.transAxes, color = 'g', fontsize=abcd_font_size)
#leg = ax3.legend();
#ax3.legend( frameon=False,loc='bottom',ncol=2,prop={'size': legend_font_size})

ax4= fig3.add_subplot(grid[2, 1])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Kamino)+1) )
for i in range(len(z0First_space_Kamino)):
    ax4.plot(FC_space_Kamino,PkPrm_Kamino[i,:],color=colors[i], linewidth=trace_width,label='Prime Conc.'+str(z0First_space_Kamino[i]))
ax4.set_ylim([0,0.3]); # ax4.set_xlim([0,10**2])      
ax4.set_xscale('log')
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('Kamino 2017',  color=mycolors[7],
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax4.text(-0.1, 1.18, 'F', ha='center',va='center',
         transform = ax4.transAxes, color = 'g', fontsize=abcd_font_size)
#leg = ax4.legend();
#ax4.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})

fig3.text(0.5, 0.04, r'$cAMP_{e}$'+' Fold Change', ha='center', va='center',fontsize=label_font_size)
fig3.text(0.05, 0.4, 'Second Spike Prominence, normalized to 1st spike height A.U.', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

#%% Plot all 4 models- 2020/3/22, with Maeda, no Sgro w/o noise
abcd_font_size = 28
title_font_size = 20
label_font_size = 18 # 16
sublabel_font_size = 18
tick_font_size = 18
legend_font_size = 14
trace_width = 3

#abcd_font_size = 28
#label_font_size=24
#title_font_size = 26
#sublabel_font_size = 22
#trace_width=3
#tick_font_size=20

mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00'] 
          
#fig3 = plt.figure(figsize=(16, 14))
#grid = plt.GridSpec(2, 2, wspace=0.15, hspace=0.25)
fig3 = plt.figure(figsize=(9,15))
grid = plt.GridSpec(3, 2, wspace=0.3, hspace=0.5)

ax0u = fig3.add_axes([0.12, 0.88, 0.305,0.075],xticks=[0,10,20,30],yticks=[0,2])
ax0u.plot(t_plot_Kamino_B,signal_trace, 'k',linewidth=trace_width)
ax0u.set_xlim([0,30]); ax0u.set_ylim([-0.5,3.5]); 
ax0u.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0u.text(-0.25 , 1.3, 'A', ha='center',va='center',
     transform = ax0u.transAxes, color = 'b', fontsize=abcd_font_size)
ax0u.text(-0.28, 0.5, r'$cAMP_{e}$'+'\n input',ha='center',va='center', rotation='vertical',
     transform = ax0u.transAxes, color = 'k', fontsize=sublabel_font_size)

ax0l = fig3.add_axes([0.12, 0.77, 0.305,0.075],xticks=[0,10,20,30],yticks=[0,0.2,0.4])
ax0l.plot(t_plot_Kamino_B, y_trace_B, linewidth=trace_width)
ax0l.set_xlim([0,30]); ax0l.set_ylim([-0.05,0.4]); 
ax0l.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0l.set_xlabel('Time, A.U.', size=sublabel_font_size)
ax0l.text(-0.28, 0.5, r'$cAMP_{i}$'+'\n response',ha='center',va='center', rotation='vertical',
     transform = ax0l.transAxes, color = 'k', fontsize=sublabel_font_size)
#ax0l.annotate('annotate', xy=(0.1,0.01), xytext=(1,1),rotation = 'vertical',
#            arrowprops={'arrowstyle': '<->'}, va='center')

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
ax0.text(-0.27, 0.5, 'Second peak \n prominence, A.U.', ha='center',va='center', rotation ='vertical',
         transform = ax0.transAxes, fontsize=label_font_size)



ax1= fig3.add_subplot(grid[1, 0])
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
ax1.set_title('Martiel 1987', color=mycolors[0],
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax1.text(-0.1, 1.15, 'B', ha='center',va='center',
         transform = ax1.transAxes, color = 'g', fontsize=abcd_font_size)
#leg = ax1.legend();
#ax1.legend( frameon=False,loc='bottom center',ncol=1,prop={'size': legend_font_size})

ax2= fig3.add_subplot(grid[1, 1])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Loomis)+1) )
for i in range(len(z0First_space_Loomis)):
    ax2.plot(FC_space_Loomis,PkPrm_Loomis[i,:],color=colors[i], linewidth=trace_width,label='Prime Conc.'+str(z0First_space_Loomis[i]))

# ax2.set_ylim([0,0.6]); 
ax2.set_xlim([FC_space_Loomis[0],FC_space_Loomis[-1]])
ax2.set_xscale('log')
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Maeda 2004', color=mycolors[1],
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax2.text(-0.1, 1.18, 'C', ha='center',va='center',
         transform = ax2.transAxes, color = 'g', fontsize=abcd_font_size)

ax4= fig3.add_subplot(grid[2, 1])
colors = plt.cm.summer(np.linspace(0,1,len(z0First_space_Kamino)+1) )
for i in range(len(z0First_space_Kamino)):
    ax4.plot(FC_space_Kamino,PkPrm_Kamino[i,:],color=colors[i], linewidth=trace_width,label='Prime Conc.'+str(z0First_space_Kamino[i]))
ax4.set_ylim([0,1.1]); 
ax4.set_xlim([FC_space_Kamino[0],FC_space_Kamino[-1]])    
ax4.set_xscale('log')
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('Kamino 2017',  color=mycolors[7],
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax4.text(-0.1, 1.18, 'D', ha='center',va='center',
         transform = ax4.transAxes, color = 'g', fontsize=abcd_font_size)
#leg = ax4.legend();
#ax4.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': legend_font_size})

# Sgro with noise
ax5= fig3.add_subplot(grid[2, 0])
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
ax5.set_title('Sgro 2015', color=mycolors[5],
              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax5.text(-0.1, 1.18, 'E', ha='center',va='center',
         transform = ax5.transAxes, color = 'g', fontsize=abcd_font_size)
#leg = ax3.legend();
#ax3.legend( frameon=False,loc='bottom',ncol=2,prop={'size': legend_font_size})

## Sgro w/o noise
#ax3= fig3.add_subplot(grid[2, 1])
#for i in range(len(z0First_space_Sgro)):
#    ax3.plot(FC_space_Sgro_noise,PkPrm_Sgro_woNoise_norm[i,:],color=colors[i], linewidth=trace_width,label='Prime Conc.'+str(z0First_space_Sgro_noise[i]))
#
#ax3.set_ylim([-0.1,1]); # ax3.set_xlim([0,10**1.5])
##ax3.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
##ax3.set_xlabel(r'$cAMP_{ext}$'+' fold change',fontsize=label_font_size)
#ax3.set_xscale('log')
#ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
#ax3.set_title('Sgro 2015 (w/o noise)',  color=mycolors[5],
#              fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
#ax3.text(-0.1, 1.18, 'E', ha='center',va='center',
#         transform = ax3.transAxes, color = 'g', fontsize=abcd_font_size)
##leg = ax3.legend();
##ax3.legend( frameon=False,loc='bottom',ncol=2,prop={'size': legend_font_size})



fig3.text(0.5, 0.04, r'$cAMP_{e}$'+' fold change', ha='center', va='center',fontsize=label_font_size)
fig3.text(0.03, 0.4, 'Second peak prominence', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
# fig3.text(0.03, 0.4, 'Second Spike Prominence, A.U.\n normalized to 1st spike,', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
#%% Save all outputs in npz file
np.savez('single_cell_FCD_20322.npz', FC_space_Gold = FC_space_Gold,z0First_space_Gold = z0First_space_Gold, PkPrm_Gold=PkPrm_Gold,PkPrm_Gold_norm=PkPrm_Gold_norm,
         FC_space_Sgro = FC_space_Sgro, z0First_space_Sgro = z0First_space_Sgro, PkPrm_Sgro_woNoise=PkPrm_Sgro_woNoise,PkPrm_Sgro_woNoise_norm=PkPrm_Sgro_woNoise_norm,
         FC_space_Sgro_noise = FC_space_Sgro_noise, z0First_space_Sgro_noise = z0First_space_Sgro_noise, 
         PkPrm_Sgro_noise=PkPrm_Sgro_noise, PkPrm_Sgro_noise_norm=PkPrm_Sgro_noise_norm,
         PkPrm_Sgro_mean_noise=PkPrm_Sgro_mean_noise, PkPrm_Sgro_se_noise = PkPrm_Sgro_se_noise,  
         FC_space_Kamino = FC_space_Kamino, z0First_space_Kamino = z0First_space_Kamino, 
         PkPrm_Kamino=PkPrm_Kamino, PkPrm_Kamino_norm=PkPrm_Kamino_norm,
         FC_space_Loomis = FC_space_Loomis, z0First_space_Loomis = z0First_space_Loomis, 
         PkPrm_Loomis=PkPrm_Loomis, PkPrm_Loomis_norm=PkPrm_Loomis_norm)

