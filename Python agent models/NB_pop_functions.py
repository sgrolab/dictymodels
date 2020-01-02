# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:23:10 2019

@author: Chuqiao Huyan

Functions needed top run the single cell (SC) iPython notebooks: 
SC_spike_and_oscillations.ipynb
SC_step_ramp_input.ipynb

"""

import numpy as np
import random
import math
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import chirp, find_peaks, peak_widths
import pandas as pd
import scipy.io

from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_pop_3var
from MaedaLoomis2004_agent_and_pop_FUN import MaedaLoomis2004_pop
from Gregor2010_agent_and_pop_FUN import Gregor2010_pop
from Sgro2015_agent_and_pop_FUN import Sgro2015_pop
from Kamino2017_agent_and_pop_FUN import Kamino2017_pop

Nt_Gregor = 6 
Nt_Sgro = 27
Nt_Goldbeter = 6.94
Nt_Maeda = 3.57
Nt_Kamino = 5.22
#%%
# Goldbeter 1986 population response simulation
def Goldbeter_pop(Goldbeter3PopParam,dt,t,campExt_influx_trace):
    # Initializations
    p0=0.8; a0=3; b0=0.9; g0=0
    Goldbeter3_pop= Goldbeter1987_pop_3var([1,1],[p0,a0,b0,g0],Goldbeter3PopParam)
    p_trace=[p0]; b_trace=[b0]; g_trace=[g0]
    for i in range(len(t)-1):
        p_next,b_next,g_next= Goldbeter3_pop.update(dt,a0,campExt_influx_trace[i])
        p_trace.append(p_next)
        b_trace.append(b_next)
        g_trace.append(g_next)
    # Convert into np array
    b_trace = np.array(b_trace);
    p_trace = np.array(p_trace); 
    t_plot_Goldbeter = np.array(t)
    return t_plot_Goldbeter, b_trace, p_trace

# Maeda 2004 population response simulation
def Maeda_pop(MaedaPopParam,dt,t,campExt_influx_trace):
    # initializations
    ACA0=0.1; PKA0=0.1; ERK20=0.1; RegA0=0.1; cAMPi0=0.01; 
    cAMPe0=0.1; CAR10=0.1
    state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]
    Maeda_agent=MaedaLoomis2004_pop([1,1],state0,MaedaPopParam)
    ACA_trace=[ACA0]; PKA_trace=[PKA0]; ERK2_trace=[ERK20]
    RegA_trace=[RegA0]; cAMPi_trace=[cAMPi0];  cAMPe_trace=[cAMPe0]
    CAR1_trace=[CAR10]
    for i in range(len(t)-1):     
        ACA_next,PKA_next,ERK2_next,RegA_next,\
        cAMPi_next,cAMPe_next,CAR1_next=Maeda_agent.update(dt,campExt_influx_trace[i])       
        ACA_trace.append(ACA_next)
        PKA_trace.append(PKA_next)
        ERK2_trace.append(ERK2_next)
        RegA_trace.append(RegA_next)
        cAMPi_trace.append(cAMPi_next)
        cAMPe_trace.append(cAMPe_next)
        CAR1_trace.append(CAR1_next)    
    ERK2_trace = np.array(ERK2_trace) # vectorize p_trace
    cAMPi_trace = np.array(cAMPi_trace)
    t_plot_Maeda = np.array(t)
    return t_plot_Maeda, cAMPi_trace, ERK2_trace

# Gregor 2010 population response simulation
def Gregor_pop(GregorPopParam,dt,t, campExt_influx_trace, time_separation = 0):
    
    # Initializations
    Amax=GregorPopParam['Amax']   
    Abas=GregorPopParam['Abas']  
    eta = GregorPopParam['eta']
    rho = Amax=GregorPopParam['rho'] 
    k = Amax=GregorPopParam['k'] 
    Vt = Amax=GregorPopParam['Vt'] 
    Nc = Amax=GregorPopParam['Nc'] 
    campCyto0 = 7.5*np.ones(Nc)
    sinthetai0 = (campCyto0*2-Amax-Abas)/(-Amax+Abas)
    thetai0 = np.arcsin(sinthetai0)
    campExt0 = 0 
    
    Gregor_pop_obj=Gregor2010_pop(campCyto0, thetai0, campExt0, GregorPopParam)
    gregor_thetai_trace=np.zeros((Nc,len(t))) ;  gregor_thetai_trace[:,0] = thetai0; 
    gregor_campCyto_trace=np.zeros((Nc,len(t))) ; gregor_campCyto_trace[:,0] = campCyto0
    gregor_campExt_trace=np.zeros(len(t)) 
    
    for i in range(len(t)-1):
        thetai_next, campCyto_next, campExt_next = Gregor_pop_obj.update(dt,eta,rho,k,Vt,time_separation,campExt_influx_trace[i])
        gregor_thetai_trace[:,i+1] = thetai_next
        gregor_campCyto_trace[:,i+1] = campCyto_next
        gregor_campExt_trace[i+1] = campExt_next
    #Traces
    gregor_thetai_trace= np.array(gregor_thetai_trace) 
    gregor_campCyto_trace= np.array(gregor_campCyto_trace) 
    t_plot_Gregor = np.array(t)
    return t_plot_Gregor,  gregor_campCyto_trace, gregor_thetai_trace

# Gregor 2010 population response simulation
def Sgro_pop(SgroPopParam,dt,t,cAMPext_influx_trace, time_separation = 0):
    N = SgroPopParam['N']
    Na = SgroPopParam['Na']
    A_trace_offset = SgroPopParam['offset_A']

    # create an object and initializations
    A0=-1.5; R0=-0.5; cAMPext0 = 0
    Sgro_pop_obj=Sgro2015_pop(A0,R0,cAMPext0, SgroPopParam)
    A_trace_orig = np.zeros(N,len(t));  A_trace_orig[:,0] = A0
    A_trace_plot = np.zeros(N,len(t));  A_trace_plot[:,0] = (A0+A_trace_offset)/Na
    R_trace_orig = np.zeros(N,len(t));  A_trace_orig[:,0] = A0
    cAMPext_trace = np.zeros((len(t),1)); cAMPext_trace[0] = cAMPext0
    for i in range(len(t)-1):
        A_next,R_next, cAMPext_next = Sgro_pop.update(dt, time_separation,cAMPext_influx_trace[i])
        A_trace_orig[:,i] = A_next
        R_trace_orig[:,i] = R_next
        cAMPext_trace[i] = cAMPext_next
        
    # Traces
    A_trace_offset=1.5
    A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
    A_trace_plot=(A_trace_orig+A_trace_offset)/SgroAgentParam['Na']
    R_trace_orig = np.array(R_trace_orig)
    t_plot_Sgro = np.array(t)
    return t_plot_Sgro, A_trace_plot,  R_trace_orig

def plot_POP_oscillation(t_plot,pop_trace_plot, cAMPext_influx, t_tot, stim_time, SC_traces = 0, SC_traces_idx = 0):
    fig = plt.figure(figsize=(11, 3))
    grid = plt.GridSpec(1, 1, wspace=0.5, hspace=0.3)
    ax1= fig.add_subplot(grid[0, 0])
    # plot defined single cell traces
    if SC_traces_idx != 0:
        for this_idx in SC_traces_idx:
            this_trace = SC_traces[this_idx,:]/np.amax(SC_traces[this_idx,:])
            ax1.plot(t_plot,this_trace, color='k',alpha=0.6, linewidth=2)     
    # Plot population mean
    ax1.plot(t_plot, pop_trace_plot, color='g' ,linewidth=3)
    ax1.set_title(r'$cAMP_{e}$ input = '+str(cAMPext_influx), fontsize=20)
    ax1.set_xlabel('Time, A.U.', fontsize=20); 
    ax1.set_ylabel(r'$cAMP_{i}$ response',fontsize=20 )
    ax1.tick_params(grid_linewidth = 15, labelsize = 20)
    ax1.axvspan(t_tot*stim_time, t_tot, alpha=0.2, color='g')
#    ax1.set_xlim([0,30]); ax1.set_ylim([-0.25,1.25])
    plt.show()




def Kamino2017_SC(KaminoAgentParam,dt,t,signal_trace):
    # Initializations
    x0=0.01; y0=0.05; z0=0.005
    Kamino_agent=Kamino2017_agent([x0,y0,z0],KaminoAgentParam)
    x_trace=[x0]; y_trace=[y0] 
    for i in range(len(t)-1):
        x_next,y_next,z_next= Kamino_agent.update( dt, signal_trace[i])
        x_trace.append(x_next)
        y_trace.append(y_next)               
    # Convert into np array
    x_trace = np.array(x_trace) 
    y_trace = np.array(y_trace)
    t_plot_Kamino = np.array(t)
    return t_plot_Kamino, y_trace, x_trace


def SC_FCD (z0First_space, FC_space, cAMP, Nt, dt, t,
            prm_lims, SecondPkEndTime, single_trace_to_plot,
            SC_model, AgentParam):
    # Initialize peak prominence array
    PkPrm = np.zeros((len(z0First_space), len(FC_space))) 
    
    for j in range(len(z0First_space)):
        z0First = z0First_space[j]
        signal_trace=z0First*np.ones(len(t))
        for k in range(len(FC_space)):
            FC = FC_space[k]
            stim_time_step=int(round(0.5*max(t)/dt)) # at this time second step input is applied
            signal_trace[stim_time_step:] = FC * z0First
            t_plot, cAMPi_trace,__ = SC_model( AgentParam,dt, t, signal_trace)
            
            # cAMPi_trace_second=cAMPi_trace[stim_time_step:] # the second part of trace, second spike
            # end_time_step= stim_time_step + int(round(2.5*Nt/dt))
            end_time_step= int(math.floor( SecondPkEndTime *Nt/dt))
            cAMPi_trace_second=cAMPi_trace[stim_time_step:end_time_step]
            PkPos, PkProperties = find_peaks(cAMPi_trace_second, prominence=(prm_lims[0],prm_lims[1]))
            # Check find_peaks
            # plt.plot(cAMPi_trace_second)
            # plt.plot(peaks, cAMPi_trace_second[peaks], "x")
            if PkPos.size: # if there is a second spike
                PkPrm[j,k]=PkProperties["prominences"][0]
            else:
                PkPrm[j,k]=0 # if there is no second spike
            
                        # plot single cell traces as selected 
            if j in single_trace_to_plot[:, 0] and k in single_trace_to_plot[:, 1]:
                            fig = plt.figure(figsize=(5,3)); grid = plt.GridSpec(3, 1,hspace= 0.3)
                            ax1= fig.add_subplot(grid[0, 0])
                            ax1.plot(t_plot,signal_trace)
                            ax1.set_ylabel(r'$cAMP_{e}$'); 
                            ax2= fig.add_subplot(grid[1:, 0])
                            ax2.plot(t_plot,cAMPi_trace)
                            ax2.set_ylabel(r'$cAMP_{i}$'); ax2.set_xlabel('Time, A.U.');
#                            ax2.set_title('Priming conc.  '+str(z0First)+' fold change '+ str(FC)+
#                                          '\n 2nd peak prominence='+str(PkPrm[j,k]))
                            ax1.set_title('Priming conc.  '+'{:#.2n}'.format(np.float64(z0First)) +
                                ', fold change '+ '{:#.2n}'.format(np.float64(FC)) + 
                                '\n 2nd peak prominence='+ '{:#.2n}'.format(np.float64(PkPrm[j,k])))
                            plt.show()
    return PkPrm

def SC_FCD_plot(z0First_space, FC_space, PkPrm):
    # plot FC vs. second response amplitude
    colors = plt.cm.summer(np.linspace(0,1,len(z0First_space)+1))    
    fig3 = plt.figure(figsize=(6, 6))
    grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
    
    ax1= fig3.add_subplot(grid[0, 0])
    
    for i in range(len(z0First_space)):
        ax1.plot(FC_space,PkPrm[i,:], color = colors[i], label='Priming cAMP='+str(z0First_space[i]))
    
    ax1.set_ylabel( 'Second spike prominence')
    ax1.set_xlabel('Fold change in extracellular cAMP')
    ax1.set_xscale('log')
    ax1.set_title('Second spike prominence')
    leg = ax1.legend();  

def SC_FCD_Noise (z0First_space,  FC_space,num_of_runs, cAMP, Nt, dt, t,
            prm_lims, SecondPkEndTime, single_trace_to_plot,
            SC_model, AgentParam):
    # make an array of all the runs
    run_time_space =np.arange(0,num_of_runs,1)
    # Initialize 2nd peak prominence array
    PkPrm_noise = np.zeros((len(z0First_space), len(FC_space), len(run_time_space)))
    
    for j in range(len(z0First_space)):
        z0First = z0First_space[j]
        signal_trace=z0First*np.ones(len(t))
        for k in range(len(FC_space)):
            FC = FC_space[k]
            stim_time_step=int(round(0.5*max(t)/dt)) # at this time second step input is applied
            signal_trace[stim_time_step:] = FC * z0First
            # cAMPi response traces of all the test runs
            cAMPi_noise_all_traces = np.zeros([num_of_runs, len(t)])
            
            for test in run_time_space:
                    
                t_plot, cAMPi_noise,__ = SC_model( AgentParam,dt, t, signal_trace)
                cAMPi_noise_all_traces[test,:] = cAMPi_noise
                end_time_step= int(math.floor( SecondPkEndTime *Nt/dt))
                cAMPi_trace_second=cAMPi_noise[stim_time_step:end_time_step]
                PkPos, PkProperties = find_peaks(cAMPi_trace_second, prominence=(prm_lims[0],prm_lims[1]))
                # Check find_peaks
                # plt.plot(cAMPi_trace_second)
                # plt.plot(peaks, cAMPi_trace_second[peaks], "x")
                if PkPos.size: # if there is a second spike
                    PkPrm_noise[j,k,test]=PkProperties["prominences"][0]
                else:
                    PkPrm_noise[j,k,test]=0 # if there is no second spike
            # optional: plot cAMPi response traces of all the runs    
            if j in single_trace_to_plot[:, 0] and k in single_trace_to_plot[:, 1]:
                fig = plt.figure(figsize=(5,3)); grid = plt.GridSpec(3, 1,hspace= 0.3)
                ax1= fig.add_subplot(grid[0, 0])
                ax1.plot(t_plot,signal_trace)
                ax1.set_ylabel(r'$cAMP_{e}$')
                ax1.set_title('Priming conc.  '+'{:#.2n}'.format(np.float64(z0First)) +
                    ', fold change '+ '{:#.2n}'.format(np.float64(FC)))   
                ax2= fig.add_subplot(grid[1:, 0])
                ax2.set_ylabel(r'$cAMP_{i}$'); ax2.set_xlabel('Time, A.U.');
                for test in run_time_space:
                    ax2.plot(t_plot,cAMPi_noise_all_traces[test,:])
                plt.show()
            
        print('The '+str(j)+'th priming concentration is finished')
    
    PkPrm_mean_noise=np.mean(PkPrm_noise,axis=2)   
    PkPrm_se_noise = np.std(PkPrm_noise,axis=2)/math.sqrt(num_of_runs)
    # PkPrm_std_noise = np.std(PkPrm_noise,axis=2)
    
        # plot single cell traces as selected 
            
    return PkPrm_noise, PkPrm_mean_noise, PkPrm_se_noise 

def SC_FCD_Noise_plot(z0First_space, FC_space, PkPrm_mean_noise,PkPrm_se_noise):
    colors = plt.cm.summer(np.linspace(0,1,len(z0First_space)+1))
    fig3 = plt.figure(figsize=(8, 7))
    grid = plt.GridSpec(1, 1, wspace=0.1, hspace=0.35)
    ax3= fig3.add_subplot(grid[0, 0])
    for i in range(len(z0First_space)):
        ax3.plot(FC_space, PkPrm_mean_noise[i,:],'o-', color = colors[i], lw = 3, ms = 10,  label='Priming cAMP='+str(z0First_space[i]))
        ax3.errorbar(FC_space, PkPrm_mean_noise[i,:], yerr=PkPrm_se_noise[i,:],
                     fmt = 'o', color=colors[i], ecolor= colors[i], elinewidth=3, capsize=10, capthick=3)
    ax3.set_ylim([-0.1,1])
    ax3.set_ylabel( '2nd spike prominence',fontsize=20)
    ax3.set_xlabel(r'$cAMP_{e}$'+' fold change',fontsize=20)
    ax3.set_xscale('log')
    ax3.tick_params(axis='both', which='major', labelsize=16)
    ax3.set_title('2d peak prominence, with noise', fontdict={'fontsize': 20, 'fontweight': 'medium'})
    leg = ax3.legend();
    ax3.legend( frameon=False,loc='bottom center',ncol=2,prop={'size': 15})
    



def SC_entrainment(period_space, PkWdth_space, NumofCycle, cAMP, Nt, dt,  
                   single_trace_to_plot,
                   SC_model, AgentParam):
    # Initialize 
    MeanR = np.empty((len(PkWdth_space), len(period_space))) 
    MeanR[:] = np.nan # Mean cross-correlation of spike trace in cycle of stimulus, as nan matrix
    
    for  j in range(len(period_space)):# period in period_space
        period = period_space[j]
        for k in range(len(PkWdth_space)): # PkWdth in PkWdth_space
            PkWdth = PkWdth_space[k]
            if PkWdth>=period: 
                pass
            else:
                stim_1cycle=np.zeros(int(period*Nt/dt))
                stim_1cycle[0:int(PkWdth*Nt/dt)]=cAMP;
                signal_trace=np.concatenate((np.zeros(int(1*Nt/dt)), np.tile(stim_1cycle,(NumofCycle)),np.zeros(int(1*Nt/dt))),axis=0);
                
                t = list(np.arange(0,len(signal_trace)*dt,dt))
                
                t_plot, cAMPi_trace,__ = SC_model( AgentParam,dt, t, signal_trace)
                t_plot=np.array(t)/Nt
                r=np.zeros(NumofCycle-1) # list that stores correlation coefficient to the first peak
                InitPk=cAMPi_trace[int(1*Nt/dt): int((1+period)*Nt/dt)]
                for m in range(NumofCycle-1):
                    FollowPk=cAMPi_trace[int((1+(m+1)*period)*Nt/dt) : (int((1+(m+1)*period)*Nt/dt)+len(InitPk))]
                    R = np.corrcoef(InitPk,FollowPk)
                    r[m]=R[1,0]                        
                MeanR[k,j] = np.mean(r)
                if len(single_trace_to_plot) != 0:
                    if j in single_trace_to_plot[:, 0] and k in single_trace_to_plot[:, 1]:
                        # check signal_trace & response trace
                        fig = plt.figure(figsize=(4,3)); grid = plt.GridSpec(1, 1,hspace= 0.3)
#                        grid = plt.GridSpec(3, 1,hspace= 0.3)
#                        ax1= fig.add_subplot(grid[0, 0])
#                        ax1.plot(t_plot,signal_trace)
#                        ax1.set_ylabel(r'$cAMP_{e}$'); 
                        ax2= fig.add_subplot(grid[0, 0])
                        ax2.plot(t_plot,cAMPi_trace)
                        for i in range(NumofCycle):
                            ax2.axvspan(1+i*period, 1+PkWdth + i*period, alpha=0.2, color='g')
                        ax2.set_ylabel(r'$cAMP_{i}$'); ax2.set_xlabel('Time, A.U.');
                        ax2.set_title('period '+'{:#.2n}'.format(period)+' stim peak width '+ '{:#.2n}'.format(PkWdth)+
                            '\n Entrainment quality='+ '{:#.2n}'.format(MeanR[k,j]))
                        plt.show()
    return MeanR

def SC_entrainment_plot(period_space, PkWdth_space,MeanR,clim):                  
    #% plot period- peak width- entrainment quality heatmap
    label_font_size = 16
    tick_font_size = 16
    # fig3 = plt.figure(figsize=(7, 6))
    fig3 = plt.figure(figsize=(5, 5))
    grid = plt.GridSpec(1, 1, wspace=0.1, hspace=0.5)    
    ax1= fig3.add_subplot(grid[0, 0])
    heatmap = ax1.pcolor(period_space, PkWdth_space,MeanR, cmap='jet')
    heatmap.set_clim(clim[0],clim[1])
    fig3.colorbar(heatmap, ax=ax1)
    ax1.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
    ax1.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
    ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
    plt.show()


def plot_SC_spike_oscillation(t_plot,signal_trace,trace_plot,step_signal, response_name = r'$cAMP_{i}$'):
    # Initialize     
    fig5 = plt.figure(figsize=(6, 6))
    grid = plt.GridSpec(3, 1)
    
    ax1= fig5.add_subplot(grid[0, 0])
    ax1.plot(t_plot,signal_trace)
    ax1.set_ylabel( r'$cAMP_{e}$ input' )
    # ax1.set_xlabel('Time',fontsize=label_font_size)
    ax1.set_title(r'$cAMP_e$ input level '+str(step_signal))
    
    ax2= fig5.add_subplot(grid[1:, 0])
    ax2.plot(t_plot,trace_plot, color='g')
    ax2.set_ylabel(response_name)
    # ax2.yaxis.label.set_color('g')
    ax2.set_xlabel('Time')
    plt.show()

def plot_SC_step_ramp(t_plot, input_trace, response_trace,response_name =  r'$cAMP_{i}, A.U.$'):
    
    fig = plt.figure(figsize=(6,6))
    grid = plt.GridSpec(3, 1)
    ax0 = fig.add_subplot(grid[0, 0])
    ax0.plot(t_plot,input_trace,color='k')
    ax0.set_ylabel(r'$cAMP_{e}$(nM)')
    ax0.set_xlabel('Time, A.U.')
    
    ax1= fig.add_subplot(grid[1:, 0])
    ax1.plot(t_plot,response_trace,color='darkgreen')
    # ax1.axvspan(2.5, 7.5, alpha=0.2, color='b'); ax1.axvspan(12.5, 17.5, alpha=0.2, color='g')
    ax1.set_ylabel(response_name)
    ax1.set_xlabel('Time (min)')
    ax1.tick_params(axis='both', which='major')
    plt.show()

