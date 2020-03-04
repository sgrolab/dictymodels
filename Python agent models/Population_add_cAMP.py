# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 2019

@author: Chuqiao

Populationoscillations with added [cAMP]ext

"""
import pandas as pd
import numpy as np
import random
import math
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import find_peaks
from time import perf_counter 
import scipy.io

# Normalization parameters
from NormParam import *
#%% Experimental data
my_dir = r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/'

Sgro2015Figure6excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheetname='Figure6')

mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']     
          
# label_font_size = 18; trace_width = 3; tick_font_size=15

abcd_font_size = 28
label_font_size=24
title_font_size = 26
sublabel_font_size = 22
trace_width=3
tick_font_size=20

fig5 = plt.figure(figsize=(11, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.3)

#void384 data = MATLAB structure, need to figure out what is what

#ax01= fig5.add_subplot(grid[0, 0])
#ax01.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["No External cAMP Mean Trace"],
#                               color = 'k', linewidth=trace_width)
#ax01.axvspan(60, 120, alpha=0.2, color='b')
#ax01.set_ylim([-0.1,0.6]);ax01.set_xlim([0,120])
#ax01.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
#ax01.text(0.7,0.85,' No External cAMP', horizontalalignment='center',verticalalignment='center',
#         transform = ax01.transAxes, color = 'k', fontsize=tick_font_size)

ax02= fig5.add_subplot(grid[0, 0])
ax02.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["Low External cAMP Mean Trace"],
                               color = 'k', linewidth=trace_width)
ax02.axvspan(60, 120, alpha=0.2, color='b')
ax02.set_ylim([-0.1,0.6]);ax02.set_xlim([0,120])
ax02.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax02.text(0.7,0.75,' Low External cAMP, \n 5-10nM', horizontalalignment='center',verticalalignment='center',
         transform = ax02.transAxes, color = 'k', fontsize=tick_font_size)
ax02.set_title('Experiment',size=title_font_size)
ax03= fig5.add_subplot(grid[1, 0])
ax03.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["Intermediate External cAMP Mean Trace"],
                               color = 'k', linewidth=trace_width)
ax03.axvspan(60, 120, alpha=0.2, color='b')
ax03.set_ylim([-0.1,0.6]);ax03.set_xlim([0,120])
ax03.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax03.text(0.7,0.75,' Intermediate External \n cAMP, 10-20nM', horizontalalignment='center',verticalalignment='center',
         transform = ax03.transAxes, color = 'k', fontsize=tick_font_size)

ax04= fig5.add_subplot(grid[2, 0])
ax04.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["High External cAMP Mean Trace"],
                               color = 'k', linewidth=trace_width)
ax04.axvspan(60, 120, alpha=0.2, color='b')
ax04.set_ylim([-0.1,0.6]);ax04.set_xlim([0,120])
ax04.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax04.text(0.7,0.75,' High External cAMP, \n 100nM', horizontalalignment='center',verticalalignment='center',
         transform = ax04.transAxes, color = 'k', fontsize=tick_font_size)

fig5.text(0.02, 0.9, 'A', color='b', fontsize=abcd_font_size, ha='center')

fig5.text(0.02, 0.5, r'FRET Signal, A.U.', fontsize=label_font_size,va='center', rotation='vertical')
fig5.text(0.5, 0.05, 'Time (min)', fontsize=label_font_size, ha='center')

#%% Sgro 2015
from Sgro2015_agent_and_pop_FUN import Sgro2015_pop

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5; sigma = 0.15 # noise strength
N = 100 # number of cells in the population
rho = 10**(-3.5); j = 0.5
SgroPopParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':sigma,'N':N,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'flux_thrs':0, 'rho': rho,'j': j}

dt=0.005 ; t_tot=30*Nt_Sgro; t=list(np.arange(0,t_tot,dt))

# define extracellular stim trace
alphafval_arr = np.array([10,20,100]) #np.logspace(0, 2.5, num=5)  # np.array([100]) #
A_traces = np.zeros((len(alphafval_arr),len(t)))

A_traces_single_cell = np.zeros((len(alphafval_arr),N,len(t)))

time_separation = 0
count = 0
for alphafval in alphafval_arr:
    stim_time_step=int(round(0.5*t_tot/dt)) # at this time step input is applied
    alphafval_trace=np.zeros(len(t))
    alphafval_trace[stim_time_step:] = alphafval
    # initializations
    A_trace_orig=np.zeros((N,len(t)))
    R_trace_orig=np.zeros((N,len(t)))
    cAMPext_trace = np.zeros((len(t),1))
    A0=-1.5*np.ones(N) # + np.random.uniform(-sigma,sigma,N); ###########
    R0=-0.5*np.ones(N) # + np.random.uniform(-sigma,sigma,N); 
    cAMPext0 = 0
    Sgro_pop=Sgro2015_pop(A0,R0,cAMPext0, SgroPopParam)
    A_trace_orig[:,0]=A0
    R_trace_orig[:,0]=R0
    cAMPext_trace[0] = cAMPext0

    for i in range(len(t)-1):
        A_now=A_trace_orig[:,i]
        R_now=R_trace_orig[:,i]
        cAMPext_now = cAMPext_trace[i]
        
        A_next,R_next, cAMPext_next = Sgro_pop.update(dt, time_separation,alphafval_trace[i])
        A_trace_orig[:,i] = A_next
        R_trace_orig[:,i] = R_next
        cAMPext_trace[i] = cAMPext_next
    
    # Traces
    A_trace_offset=1.5
    A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
    A_trace_plot=(A_trace_orig+A_trace_offset)/Nh_Sgro;
    A_trace_mean_plot = np.mean(A_trace_plot,axis = 0)
    t_plot_Sgro = np.array(t)/Nt_Sgro
    
    A_traces_single_cell[count,:,:] = A_trace_plot
    A_traces[count,:] = A_trace_mean_plot
    count = count+1
    
#    #  check simulation traces
#    label_font_size=25; trace_width=3; tick_font_size=18
#    fig,ax = plt.subplots()
#    ax.plot(t_plot_Sgro,A_trace_mean_plot,linewidth=trace_width, label= r'activator, $cAMP_{cyto}$')
#    # ax.plot(t_plot_Sgro,A_trace_plot[2,:],linewidth=trace_width, label= r'activator, $cAMP_{cyto}$')
#    # ax.set_ylim([-0.2,1.3])
#    ax.set_xlabel('Time')
#    ax.set_ylabel('Activator')
#    ax.set_title(r'Sgro 2015 group oscillation, $alpha_f$= '+str(alphafval)+', time separation= '+str(time_separation))
#    leg = ax.legend()
#    ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
#    plt.show()
    
#%% Plot  3 traces: low, medium and high [cAMP]ext 

fig = plt.figure(figsize=(11, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.3)

ax1= fig.add_subplot(grid[0, 0])
for count in range(5):
    ax1.plot(t_plot_Sgro,A_traces_single_cell[0,count,:],
             color=mycolors[5],alpha=0.3, linewidth=trace_width-1)
ax1.plot(t_plot_Sgro,A_traces[0,:], color=mycolors[5],linewidth=trace_width)
ax1.text(0.7,0.9,r'Low $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax1.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax1.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax1.axvspan(15, 30, alpha=0.2, color='g')
ax1.set_xlim([0,30]); ax1.set_ylim([-0.25,1.25])

ax2= fig.add_subplot(grid[1, 0])
for count in range(5):
    ax2.plot(t_plot_Sgro,A_traces_single_cell[1,count,:],
             color=mycolors[5],alpha=0.3, linewidth=trace_width-1)
ax2.plot(t_plot_Sgro,A_traces[1,:], color=mycolors[5],linewidth=trace_width)
    

ax2.text(0.7,0.9,r'Intermediate $cAMP_{e}$'+' input', horizontalalignment='center',verticalalignment='center',
     transform = ax2.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax2.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax2.axvspan(15, 30, alpha=0.2, color='g')
ax2.set_xlim([0,30]); ax2.set_ylim([-0.25,1.25])

ax3= fig.add_subplot(grid[2, 0])
for count in range(5):
    ax3.plot(t_plot_Sgro,A_traces_single_cell[2,count,:],
             color=mycolors[5],alpha=0.3, linewidth=trace_width-1)
ax3.plot(t_plot_Sgro,A_traces[2,:], color=mycolors[5],linewidth=trace_width)
ax3.text(0.7,0.9,r'High $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax3.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax3.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax3.axvspan(15, 30, alpha=0.2, color='g')
ax3.set_xlim([0,30]); ax3.set_ylim([-0.25,1.25])

fig.text(0.02, 0.9, 'E', color='g', fontsize=abcd_font_size, ha='center')
fig.text(0.5, 0.04, 'Time, A.U.',fontsize=label_font_size, ha='center')
fig.text(0.02, 0.5, r'$cAMP_{i}$',fontsize=label_font_size, va='center', rotation='vertical')
fig.text(0.5, 0.9, 'Sgro 2015',color = mycolors[5],fontsize=title_font_size, ha='center')
plt.show()
#for i in range(len(alphafval_arr)):
#    ax= fig.add_subplot(grid[i, 0])
#    ax.plot(t_plot_Sgro,A_traces[i,:], color=mycolors[5],linewidth=trace_width)
#    ax.text(0.45,0.8-2.2*i,r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', horizontalalignment='center',verticalalignment='center',
#         transform = ax01.transAxes, color = 'k', fontsize=tick_font_size)
##    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
#    ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
#    ax.axvspan(15, 30, alpha=0.2, color='g')
#    ax.set_xlim([0,30])



#%% Gregor 2010

from Gregor2010_agent_and_pop_FUN import Gregor2010_pop

Amax=20;  Abas=0.4 # uM
w=2*math.pi/6 # min-1
Vc=1.1e-9 # ml
St=1.33 # cm2
Sc=1.3e-6 # cm2
K=0.0004 # uM, 400 pM
c_sec= 3.6 # min-1
c_excite=1.01 # min-1
Nc=100 # Num of cells
rho = 1/12 #1/ml
Vt = 1 #chamber size ml
k = 5 #ml/min
GregorPopParam={'Amax':Amax,'Abas':Abas,'w':w,'Vc':Vc,'St':St,'Sc':Sc,'K':K,\
            'c_sec':c_sec,'c_excite':c_excite,'Nc':Nc, 'rho':rho, 'Vt':Vt,'k':k}

dt=0.005; t_tot=30*Nt_Gregor; t=list(np.arange(0,t_tot,dt))
eta=0.002 # noise stength

ext_input_arr = np.array([0.0001,10,1000]) #np.logspace(0,3,num=5)
campCyto_traces = np.zeros((len(ext_input_arr),len(t)))
campCyto_traces_single_cell = np.zeros((len(ext_input_arr),Nc,len(t)))
time_separation = 0
count = 0
for ext_input in ext_input_arr:
    stim_time_step=int(round(0.1*t_tot/dt)) # at this time step input is applied
    ext_input_trace=np.zeros(len(t))
    ext_input_trace[stim_time_step:] = ext_input
    # Initializations
    campCyto0 = 7.5*np.ones(Nc)
    sinthetai0 = (campCyto0*2-Amax-Abas)/(-Amax+Abas)
    thetai0 = np.arcsin(sinthetai0)
    campExt0 = 0 # Vc*St/Sc*rho/K*c_sec*1/Nc*np.sum(campCyto0);
    Gregor_pop=Gregor2010_pop(campCyto0, thetai0, campExt0, GregorPopParam)
    gregor_thetai_trace=np.zeros((Nc,len(t))) 
    gregor_campCyto_trace=np.zeros((Nc,len(t))) 
    gregor_campExt_trace=np.zeros(len(t)) 

    for i in range(len(t)-1):
        thetai_now=gregor_thetai_trace[:,i]
        campCyto_now=gregor_campCyto_trace[:,i]
        campExt_now=gregor_campExt_trace[i]
        thetai_next, campCyto_next, campExt_next = Gregor_pop.update(dt,eta,rho,k,Vt,time_separation,ext_input_trace[i])
        gregor_thetai_trace[:,i+1] = thetai_next
        gregor_campCyto_trace[:,i+1] = campCyto_next
        gregor_campExt_trace[i+1] = campExt_next
        
    #Traces
    # gregor_thetai_trace= np.array(gregor_thetai_trace) 
    # gregor_campExt_trace = np.array(gregor_campExt_trace)
    gregor_campCyto_trace= np.array(gregor_campCyto_trace)/Nh_Gregor
    gregor_campCyto_trace_mean= np.mean(gregor_campCyto_trace,axis = 0)
    gregor_campCyto_trace_mean = gregor_campCyto_trace_mean
    t_plot_Gregor = np.array(t)/Nt_Gregor
    campCyto_traces[count,:] = gregor_campCyto_trace_mean
    campCyto_traces_single_cell[count,:,:] = gregor_campCyto_trace
    count = count+1
#    #  check simulation traces
#    label_font_size=25; trace_width=3; tick_font_size=18
#    fig,ax = plt.subplots()
#    ax.plot(t_plot_Gregor,gregor_campCyto_trace_mean,linewidth=trace_width, label= r'activator, $cAMP_{cyto}$')
#    # ax.plot(t_plot_Gregor,gregor_campExt_trace,linewidth=trace_width, label= r'$cAMP_{ext}$')
#    ax.set_xlabel('Time')
#    ax.set_ylabel('Activator')
#    ax.set_title(r'Gregor 2010, cAMP input='+str(ext_input)+' T_sep= '+str(time_separation))
#    leg = ax.legend()
#    ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
#    plt.show()
    
# Plot  3 traces: low, medium and high [cAMP]ext 
fig = plt.figure(figsize=(11, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.3)

ax1= fig.add_subplot(grid[0, 0])
ax1.plot(t_plot_Gregor,campCyto_traces[0,:], alpha=0.8, color=mycolors[2],linewidth=trace_width+5)  
for count in range(10):
    campCyto_traces_single_cell[0,count,:] = campCyto_traces_single_cell[0,count,:]/np.amax(campCyto_traces_single_cell[0,count,:])
    ax1.plot(t_plot_Gregor,campCyto_traces_single_cell[0,count,:],
             color='k',alpha=0.6, linewidth=trace_width-1)
 

ax1.text(0.7,0.9,r'Low $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax1.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax1.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax1.axvspan(15, 30, alpha=0.2, color='g')
ax1.set_xlim([0,30]); ax1.set_ylim([0,1.5])

ax2= fig.add_subplot(grid[1, 0])
ax2.plot(t_plot_Gregor,campCyto_traces[1,:], alpha=0.8, color=mycolors[2],linewidth=trace_width+5)
for count in range(10):
    campCyto_traces_single_cell[1,count,:] = campCyto_traces_single_cell[1,count,:]/np.amax(campCyto_traces_single_cell[1,count,:])
    ax2.plot(t_plot_Gregor,campCyto_traces_single_cell[1,count,:],
             color='k',alpha=0.6, linewidth=trace_width-1)
ax2.text(0.7,0.9,r'Intermediate $cAMP_{e}$'+' input', horizontalalignment='center',verticalalignment='center',
     transform = ax2.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax2.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax2.axvspan(15, 30, alpha=0.2, color='g')
ax2.set_xlim([0,30]); ax2.set_ylim([0, 1.5])

ax3= fig.add_subplot(grid[2, 0])
ax3.plot(t_plot_Gregor,campCyto_traces[2,:], color=mycolors[2],alpha = 0.8, linewidth=trace_width+5)
for count in range(10):
    campCyto_traces_single_cell[2,count,:] = campCyto_traces_single_cell[2,count,:]/np.amax(campCyto_traces_single_cell[2,count,:])
    ax3.plot(t_plot_Gregor,campCyto_traces_single_cell[2,count,:],
             color='k',alpha=0.6, linewidth=trace_width-1)
ax3.text(0.7,0.9,r'High $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax3.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax3.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax3.axvspan(15, 30, alpha=0.2, color='g')
ax3.set_xlim([0,30]); ax3.set_ylim([0,1.5])

fig.text(0.02, 0.9, 'D', color='g', fontsize=abcd_font_size, ha='center')
fig.text(0.5, 0.04, 'Time, A.U.',fontsize=label_font_size, ha='center')
fig.text(0.02, 0.5, r'$cAMP_{i}$',fontsize=label_font_size, va='center', rotation='vertical')
fig.text(0.5, 0.9, 'Gregor 2010',color = mycolors[2],fontsize=label_font_size, ha='center')
plt.show()

#for i in range(len(alphafval_arr)):
#    ax= fig.add_subplot(grid[i, 0])
#    ax.plot(t_plot_Gregor,campCyto_traces[i,:], color=mycolors[2],linewidth=trace_width)
#    ax.text(0.45,0.8-2.2*i,r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', horizontalalignment='center',verticalalignment='center',
#         transform = ax01.transAxes, color = 'k', fontsize=tick_font_size)
##    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
#    ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
#    ax.axvspan(15, 30, alpha=0.2, color='g')
#    ax.set_xlim([0,30])
#    ax.set_ylim([-2,27])



#%% Golbeter 1987, Table II/ Fig 3 parameters
from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_pop_3var

k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10; L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01; theta=0.01
e=  1 # 0.108 # compared to 1
q=4000
sig= 0.6 # 0.57 # compared to 0.6
v=12; k= 4 # k prime in the paper
ki=1.7 # 0.958 # compared to 1.7 
kt=0.9

kc=5.4 # 3.58 # compared to 5.4
h=5 # ratio of intracellular to extracellular volume, rho 
Goldbeter3AgentParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}

dt=0.0005; t_tot=30*Nt_Goldbeter; t=list(np.arange(0,t_tot,dt))
camp_input_Goldbeter_arr = np.array([0.01,0.02,0.1]) #np.logspace(-2, -1, num=5)  #
b_traces = np.zeros((len(camp_input_Goldbeter_arr),len(t)))

count = 0
for camp_input in camp_input_Goldbeter_arr:
    stim_time_step=int(round(0.5*t_tot/dt)) # at this time step input is applied
    camp_input_trace=np.zeros(len(t))
    camp_input_trace[stim_time_step:] = camp_input
    # initializations
    p0=0.8; a0=3; b0=0.9; g0=0
    p_trace=[p0]; b_trace=[b0]; g_trace=[g0]
    Goldbeter3_pop = Goldbeter1987_pop_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
    for i in range(len(t)-1):
        p_next,b_next,g_next= Goldbeter3_pop.update(dt,a0,camp_input_trace[i])
        p_trace.append(p_next)
        b_trace.append(b_next)
        g_trace.append(g_next)
    # Convert into np array
    b_trace = np.array(b_trace);  b_trace = b_trace/Nh_Goldbeter
    p_trace = np.array(p_trace); 
    t_plot_Goldbeter = np.array(t)/Nt_Goldbeter
    b_traces[count,:] = b_trace
    count = count+1
    
#    #  check simulation traces
#    label_font_size=25; trace_width=3; tick_font_size=18
#    fig,ax = plt.subplots()
#    ax.plot(t_plot_Goldbeter,b_trace,linewidth=trace_width, label= r'b, $cAMP_{cyto}$')
#    # ax.plot(t_plot_Goldbeter,p_trace, linewidth=trace_width,label = r'p, $R_{act}/R_{tot}$')
#    # ax.set_ylim([-0.2,1.3])
#    ax.set_xlabel('Time')
#    ax.set_ylabel('b, p')
#    ax.set_title(r'Goldber 1987 pop oscillation, $cAMP_ext$='+str(camp_input))
#    leg = ax.legend()
#    ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
#    plt.show()

#%%  Plot  3 traces: low, medium and high [cAMP]ext 
fig = plt.figure(figsize=(11, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.3)

ax1= fig.add_subplot(grid[0, 0])
ax1.plot(t_plot_Goldbeter,b_traces[0,:], color=mycolors[0],linewidth=trace_width)
ax1.text(0.7,0.9,r'Low $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax1.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax1.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax1.axvspan(15, 30, alpha=0.2, color='g')
ax1.set_xlim([0,30]); ax1.set_ylim([-0.25,1.5])

ax2= fig.add_subplot(grid[1, 0])
ax2.plot(t_plot_Goldbeter,b_traces[1,:], color=mycolors[0],linewidth=trace_width)
ax2.text(0.7,0.8,r'Intermediate $cAMP_{e}$'+'\n input', horizontalalignment='center',verticalalignment='center',
     transform = ax2.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax2.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax2.axvspan(15, 30, alpha=0.2, color='g')
ax2.set_xlim([0,30]); ax2.set_ylim([-0.25,1.5])

ax3= fig.add_subplot(grid[2, 0])
ax3.plot(t_plot_Goldbeter,b_traces[2,:], color=mycolors[0],linewidth=trace_width)
ax3.text(0.7,0.8,r'High $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax3.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax3.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax3.axvspan(15, 30, alpha=0.2, color='g')
ax3.set_xlim([0,30]); ax3.set_ylim([-0.25,1.5])

fig.text(0.02, 0.9, 'B', color='g', fontsize=abcd_font_size, ha='center')
fig.text(0.5, 0.04, 'Time, A.U.',fontsize=label_font_size, ha='center')
fig.text(0.02, 0.5, r'$cAMP_{i}$',fontsize=label_font_size, va='center', rotation='vertical')
fig.text(0.5, 0.9, 'Martiel 1987',color = mycolors[0],fontsize=label_font_size, ha='center')
plt.show()

#ax2= fig.add_subplot(grid[1, 0])
#line1=ax2.plot(t_plot_Sgro,A_traces[0,:], color='g',linewidth=trace_width)
##ax2.set_ylabel('Activator',fontsize=label_font_size)
##ax2.yaxis.label.set_color('g')
#ax2.set_title(r'$\alpha_f$= '+str(alphafval_arr[0]))

#%% Maeda & Loomis 2004
from MaedaLoomis2004_agent_and_pop_FUN import MaedaLoomis2004_pop
# parameters from Maeda & Loomis 2004 paper
k1=2.0; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=1.0; k8=1.3; k9=0.3; k10=0.8
k11=0.7; k12=4.9; k13=23; k14=4.5
MaedaAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}
## parameters from Laub & Loomis 1998 paper
#k1=1.4; k2=0.9; k3=2.5; k4=1.5; k5=0.6
#k6=0.8; k7=2.0; k8=1.3; k9=0.7; k10=1.0
#k11=0.3; k12=3.1; k13=1.8; k14=1.5
#LaubAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
#            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
#            'k13':k13,'k14':k14}

dt=0.001; t_tot = (30+30)*Nt_Maeda; t=list(np.arange(0,t_tot,dt))
camp_input_Maeda_arr = np.array([0.005, 0.05, 0.5])   # np.logspace(-1, 2, num=6) # 
cAMPi_traces = np.zeros((len(camp_input_Maeda_arr),int(len(t)/2)))
count = 0
for camp_input in camp_input_Maeda_arr:
    stim_time_step=int(round(0.75*t_tot/dt)) # at this time step input is applied
    camp_input_trace=np.zeros(len(t))
    camp_input_trace[stim_time_step:] = camp_input
    
    ACA0=0.1; PKA0=0.1; ERK20=0.1; RegA0=0.1; cAMPi0=0.01; 
    cAMPe0=0.1; CAR10=0.1
    cAMPi_trace = [cAMPi0]
    state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]
    MaedaLoomis_pop=MaedaLoomis2004_pop([1,1],state0,MaedaAgentParam)
    
    for i in range(len(t)-1):
        ACA_next,PKA_next,ERK2_next,RegA_next,\
        cAMPi_next,cAMPe_next,CAR1_next=MaedaLoomis_pop.update(dt,camp_input_trace[i])
        cAMPi_trace.append(cAMPi_next)
    cAMPi_trace_later = np.array(cAMPi_trace[int(len(t)/2):]); 
    t_plot_Maeda = np.arange(0,t_tot/2,dt)/Nt_Maeda
    cAMPi_traces[count,:] = cAMPi_trace_later/Nh_Maeda
    count = count + 1
#    #  check simulation traces
#    label_font_size=25; trace_width=3; tick_font_size=18
#    fig,ax = plt.subplots()
#    ax.plot(t_plot_Maeda,cAMPi_trace_later,linewidth=trace_width, label= r'$cAMP_{cyto}$')
#    # ax.plot(t_plot_Goldbeter,p_trace, linewidth=trace_width,label = r'p, $R_{act}/R_{tot}$')
#    # ax.set_ylim([-0.2,1.3])
#    ax.set_xlabel('Time')
#    ax.set_ylabel(r'$cAMP_{cyto}$')
#    ax.set_title(r'Maeda 2004 pop oscillation, $cAMP_ext$='+str(camp_input))
#    leg = ax.legend()
#    ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
#    plt.show()
    
#%% Plot  3 traces: low, medium and high [cAMP]ext 
fig = plt.figure(figsize=(11, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.3)

ax1= fig.add_subplot(grid[0, 0])
ax1.plot(t_plot_Maeda,cAMPi_traces[0,:], color=mycolors[1],linewidth=trace_width)
ax1.text(0.7,0.9,r'Low $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax1.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax1.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax1.axvspan(15, 30, alpha=0.2, color='g')
ax1.set_xlim([0,30]); ax1.set_ylim([0.1,0.9])

ax2= fig.add_subplot(grid[1, 0])
ax2.plot(t_plot_Maeda, cAMPi_traces[1,:], color=mycolors[1],linewidth=trace_width)
ax2.text(0.7,0.8,r'Intermediate $cAMP_{e}$'+'\n input', horizontalalignment='center',verticalalignment='center',
     transform = ax2.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax2.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax2.axvspan(15, 30, alpha=0.2, color='g')
ax2.set_xlim([0,30]); ax2.set_ylim([0.1,0.9])

ax3= fig.add_subplot(grid[2, 0])
ax3.plot(t_plot_Maeda, cAMPi_traces[2,:], color=mycolors[1],linewidth=trace_width)
ax3.text(0.7,0.8,r'High $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax3.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax3.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax3.axvspan(15, 30, alpha=0.2, color='g')
ax3.set_xlim([0,30]); ax3.set_ylim([0.1,0.9])

fig.text(0.02, 0.9, 'C', color='g', fontsize=abcd_font_size, ha='center')
fig.text(0.5, 0.04, 'Time, A.U.',fontsize=label_font_size, ha='center')
fig.text(0.02, 0.5, r'$cAMP_{i}$',fontsize=label_font_size, va='center', rotation='vertical')
fig.text(0.5, 0.9, 'Maeda 2004',color = mycolors[1],fontsize=label_font_size, ha='center')
plt.show()

#for i in range(len(camp_input_Maeda_arr)):
#    ax= fig.add_subplot(grid[i+1, 0])
#    ax.plot(t_plot_Maeda,cAMPi_traces[i,:], color='g',linewidth=trace_width)
#    #ax2.set_ylabel('Activator',fontsize=label_font_size)
#    #ax2.yaxis.label.set_color('g')
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(camp_input_Maeda_arr[i]), fontsize=label_font_size)
#    ax.tick_params(grid_linewidth = 15, labelsize = tick_font_size)


#%% Kamino 2017, fig 5D group oscillations
from Kamino2017_agent_and_pop_FUN import Kamino2017_pop
   
tau=1.5; n=2; K=4; kt=2; delta=0
gamma = 3; rho = 1 # cell density
Param={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
       'gamma':gamma,'rho':rho}
dt=0.001; t_tot=30 *Nt_Kamino; t=np.arange(0,t_tot,dt)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              

camp_input_Kamino_arr = np.array([0.001,0.002,0.1])# np.logspace(-3, -1, num=5) # 
y_traces = np.zeros((len(camp_input_Kamino_arr),int(len(t))))

count = 0
for camp_input in camp_input_Kamino_arr:
    stim_time_step=int(round(0.5*t_tot/dt)) # at this time step input is applied
    camp_input_trace=np.zeros(len(t))
    camp_input_trace[stim_time_step:] = camp_input
    # initializations
    x0=0.01; y0=0.08; z0=0.01
    y_trace=[y0]; # x_trace=[x0];  z_trace=[z0]
    Kamino_pop=Kamino2017_pop([x0,y0,z0],Param)
    for i in range(len(t)-1):
        x_next,y_next,z_next=Kamino_pop.update(camp_input_trace[i],dt)
        y_trace.append(y_next)
        
    y_trace = np.array(y_trace);  y_trace = y_trace/Nh_Kamino
    t_plot_Kamino = np.array(t)/Nt_Kamino
    y_traces[count,:] = y_trace
    count = count+1
#    #  check simulation traces
#    label_font_size=25; trace_width=3; tick_font_size=18
#    fig,ax = plt.subplots()
#    ax.plot(t_plot_Kamino,y_trace,linewidth=trace_width, label= r'y, $cAMP_{cyto}$')
#    ax.set_xlabel('Time')
#    ax.set_ylabel('y')
#    ax.set_title(r'Kamino 2017 pop oscillation, $cAMP_ext$='+str(camp_input))
#    leg = ax.legend()
#    ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
#    plt.show()
# Plot  3 traces: low, medium and high [cAMP]ext 

fig = plt.figure(figsize=(11, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.3)

ax1= fig.add_subplot(grid[0, 0])
ax1.plot(t_plot_Kamino, y_traces[0,:], color=mycolors[7],linewidth=trace_width)
ax1.text(0.7,0.9,r'Low $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax1.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax1.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax1.axvspan(15, 30, alpha=0.2, color='g')
ax1.set_xlim([0,30]); ax1.set_ylim([-0.2,1.2])

ax2= fig.add_subplot(grid[1, 0])
ax2.plot(t_plot_Kamino, y_traces[1,:], color=mycolors[7],linewidth=trace_width)
ax2.text(0.7,0.8,r'Intermediate $cAMP_{e}$'+'\n input', horizontalalignment='center',verticalalignment='center',
     transform = ax2.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax2.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax2.axvspan(15, 30, alpha=0.2, color='g')
ax2.set_xlim([0,30]); ax2.set_ylim([-0.2,1.2])

ax3= fig.add_subplot(grid[2, 0])
ax3.plot(t_plot_Kamino, y_traces[2,:], color=mycolors[7],linewidth=trace_width)
ax3.text(0.7,0.8,r'High $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax3.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax3.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax3.axvspan(15, 30, alpha=0.2, color='g')
ax3.set_xlim([0,30]); ax3.set_ylim([-0.2,1.2])

fig.text(0.02, 0.9, 'F', color='g', fontsize=abcd_font_size, ha='center')
fig.text(0.5, 0.04, 'Time, A.U.',fontsize=label_font_size, ha='center')
fig.text(0.02, 0.5, r'$cAMP_{i}$',fontsize=label_font_size, va='center', rotation='vertical')
fig.text(0.5, 0.9, 'Kamino 2017',color = mycolors[7],fontsize=label_font_size, ha='center')
plt.show()
    
