# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 2020

@author: Chuqiao

Population collective oscillations, with noise in single cells

"""
import pandas as pd
import numpy as np
import random
import math
import matplotlib.pyplot as plt
# set matplotlib default font
import matplotlib
font = {'family' : 'Roboto'}
matplotlib.rc('font', **font)

from scipy import signal
from scipy.signal import find_peaks
from time import perf_counter 

from NB_pop_functions import plot_POP_oscillation

# Normalization parameters
from Params import NormParams
for key,val in NormParams.items():
        exec(key + '=val')

#%% change back to default
import matplotlib.style
import matplotlib 
matplotlib.style.use('default')
#%% Golbeter 1987,  Table II/ Fig 3 parameters

from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_pop_3var_SCnoise
from Params import Goldbeter3PopParam
sigma = 10 # noise strength
N = 100 # number of cells
kc=5.4
h=5

dt=0.001; t_tot=20; t=list(np.arange(0,t_tot*Nt_Goldbeter,dt))
nSteps = len(t)
cAMPext_influx = 0; cAMPext_influx_trace = np.ones(nSteps)*cAMPext_influx

# Initializations
p0=0.8*np.ones(N); a0=3; b0=0.9*np.ones(N); g0=0
p_trace=np.zeros((N,nSteps));p_trace[:,0] = p0; 
b_trace=np.zeros((N,nSteps));b_trace[:,0] = b0; 
g_trace=np.zeros((nSteps,1));g_trace[0] = g0; 

Goldbeter3_pop= Goldbeter1987_pop_3var_SCnoise(0,p0,b0,g0,Goldbeter3PopParam)

for i in range(len(t)-1):
    p_next,b_next,g_next= Goldbeter3_pop.update(dt,a0,cAMPext_influx_trace[i]) # cAMPext_influx_trace[i])
    p_trace[:,i+1]= p_next
    b_trace[:,i+1]= b_next
    g_trace[i+1] = g_next

# normalize
b_trace_norm = b_trace/Nh_Goldbeter
b_trace_norm_mean = np.mean(b_trace_norm,axis = 0)
t_plot_Goldbeter = np.array(t)/Nt_Goldbeter

later_portion = 0.25 # start count peaks after this X total simulation time
b_trace_norm_later = b_trace_norm[:,math.floor(nSteps * later_portion):]
b_trace_norm_mean_later=b_trace_norm_mean[math.floor(nSteps * later_portion):] # the later part of trace

t_plot_Goldbeter_short = t_plot_Goldbeter[0:math.floor(nSteps * (1-later_portion))] 

pk_find_thresh = 5/Nh_Goldbeter
pop_max = np.amax(b_trace_norm_mean_later); pop_min = np.amin(b_trace_norm_mean_later); 
PkPos, PkProperties = find_peaks(b_trace_norm_mean_later, prominence=(pk_find_thresh,pop_max))
#PkPos, PkProperties = find_peaks(b_trace_norm_mean_later, prominence=((pop_max-pop_min)*0.4,pop_max))

if len(PkPos) == 0:
        firing_rate = 0; height = 0
else: 
    firing_rate = len(PkPos)/(t_tot*(1-later_portion))

# Check find_peaks
title = 'Goldbeter1987, noise=' +str(sigma)+', density='+'{:#.3n}'.format(np.float64(1/h))+\
    ', J='+ '{:#.3n}'.format(np.float64(kc)) +'\n ,pk_find_thresh='+str(pk_find_thresh)
fig = plt.figure()
plt.plot(b_trace_norm_mean_later)
plt.plot(PkPos, b_trace_norm_mean_later[PkPos], "x")
plt.title(title+' FR= '+str(firing_rate))
    

# plot population mean and selected single cell traces
SC_traces_idx = [0,2,4,6,8,10]
plot_POP_oscillation(t_plot_Goldbeter_short,b_trace_norm_mean_later,cAMPext_influx,
                     t_tot*(1-later_portion), 0,title, b_trace_norm_later, SC_traces_idx)

##  check simulation traces
#label_font_size=25; trace_width=3; tick_font_size=18
#
#fig,ax = plt.subplots()
#ax.plot(t_plot_Goldbeter,b_trace_norm_mean,linewidth=trace_width, label= r'b, $cAMP_{cyto}$')
## ax.plot(t_plot_Goldbeter,p_trace, linewidth=trace_width,label = r'p, $R_{act}/R_{tot}$')
## ax.set_ylim([-0.2,1.3])
#ax.set_xlabel('Time')
#ax.set_ylabel('b, p')
#ax.set_title('Goldber 1987 group oscillation, noise= ')
#leg = ax.legend()
#ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
#plt.show()

## Get the oscillation period
#b_trace=np.array(b_trace) # convert list to array
#later_portion = 0.2 # start count peaks after this X total simulation time
#b_trace_later=b_trace[math.floor(len(t)*later_portion):] # the later part of trace
#PkPos, PkProperties = find_peaks(b_trace_later, prominence=(-1,2))
#
## Check find_peaks
#plt.plot(b_trace_later)
#plt.plot(PkPos, b_trace_later[PkPos], "x")
#
#Goldbeter_pop_osc_period = (1-later_portion)*t_tot / len(PkPos)
#print('group oscillation period for Golbdeter is '+str(Goldbeter_pop_osc_period))

#%% Maeda & Loomis 1998
from MaedaLoomis2004_agent_and_pop_FUN import MaedaLoomis2004_pop_SCnoise

from Params import MaedaPopParam
sigma = 0.1 # noise strength
N = 100  # number of cells
gamma = 0# [cAMP]e flow rate 
rho = 1 # cell density


dt=0.0001; t_tot=60; t=list(np.arange(0,t_tot*Nt_Maeda,dt))
nSteps = len(t)

# Initializations
ACA0=0.1*np.ones(N); PKA0=0.1*np.ones(N); ERK20=0.1*np.ones(N); 
RegA0=0.1*np.ones(N); cAMPi0=0.01*np.ones(N); cAMPe0=0.1; CAR10=0.1*np.ones(N)

ACA_trace=np.zeros((N,nSteps)); ACA_trace[:,0] = ACA0
PKA_trace=np.zeros((N,nSteps)); PKA_trace[:,0] = PKA0
ERK2_trace=np.zeros((N,nSteps)); ERK2_trace[:,0] = ERK20
RegA_trace= np.zeros((N,nSteps)); RegA_trace[:,0] = RegA0
cAMPi_trace= np.zeros((N,nSteps)); cAMPi_trace[:,0] = cAMPi0
cAMPe_trace=np.zeros((nSteps,1));cAMPe_trace[0] = cAMPe0;
CAR1_trace= np.zeros((N,nSteps)); CAR1_trace[:,0] = CAR10

MaedaLoomis_pop=MaedaLoomis2004_pop_SCnoise([1,1],ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10, MaedaPopParam)

cAMPext_influx = 0; cAMPext_influx_trace = np.ones(nSteps)*cAMPext_influx

for i in range(len(t)-1):
    ACA_next,PKA_next,ERK2_next,RegA_next,\
    cAMPi_next,cAMPe_next,CAR1_next=MaedaLoomis_pop.update(dt,cAMPext_influx_trace[i],rho,gamma)
    ACA_trace[:,i+1]= ACA_next
    PKA_trace[:,i+1]= PKA_next
    ERK2_trace[:,i+1]= ERK2_next
    RegA_trace[:,i+1]= RegA_next
    cAMPi_trace[:,i+1]= cAMPi_next
    cAMPe_trace[i+1]= cAMPe_next
    CAR1_trace[:,i+1]= CAR1_next
    

cAMPi_trace_norm = cAMPi_trace/Nh_Maeda   
cAMPi_trace_norm_mean = np.mean(cAMPi_trace_norm,axis = 0) # population mean

later_portion = 0.75 # start count peaks after this X total simulation time
cAMPi_trace_norm_later = cAMPi_trace_norm[:,math.floor(nSteps * later_portion):]
cAMPi_trace_norm_mean_later=cAMPi_trace_norm_mean[math.floor(nSteps * later_portion):] # the later part of trace

t_plot_Maeda = np.array(t)/Nt_Maeda
t_plot_Maeda_short = t_plot_Maeda[0:math.floor(nSteps * (1-later_portion))] 

pop_max = np.amax(cAMPi_trace_norm_mean_later); pop_min = np.amin(cAMPi_trace_norm_mean_later); 
pk_find_prm = 0.7
PkPos, PkProperties = find_peaks(cAMPi_trace_norm_mean_later, prominence=((pop_max-pop_min)*pk_find_prm,pop_max))

if len(PkPos) == 0:
        firing_rate = 0; height = 0
else: 
    firing_rate = len(PkPos)/(t_tot*(1-later_portion))
        
# Check find_peaks
title = 'Maeda 2004, noise=' +str(sigma)+', density='+str(rho)+', J='+str(gamma)+',pk_find_prm='+str(pk_find_prm)
fig = plt.figure()
plt.plot(cAMPi_trace_norm_mean_later)
plt.plot(PkPos, cAMPi_trace_norm_mean_later[PkPos], "x")
plt.title(title+'\n FR= '+str(firing_rate))

# plot population mean and selected single cell traces
SC_traces_idx = [0,2,4,6,8,10]
title = r'Maeda2004, $cAMP_{e}$ input = '+str(cAMPext_influx)
plot_POP_oscillation(t_plot_Maeda_short,cAMPi_trace_norm_mean_later,cAMPext_influx,
                     t_tot*(1-later_portion), 0,title, cAMPi_trace_norm_later, SC_traces_idx)


#%% Kamino 2017, fig 5D group oscillations
from Kamino2017_agent_and_pop_FUN import Kamino2017_pop_SCnoise

from Params import KaminoPopParam
KaminoPopParam['delta'] = 0 # only for population oscillations and population add cAMP

gamma = 3;  # dilution rate
rho = 1 # cell density

dt=0.001
t_tot= 20                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
t=np.arange(0,t_tot*Nt_Kamino,dt)
nSteps = len(t)
cAMPext_influx = 0
cAMPext_influx_trace = cAMPext_influx*np.ones(len(t)) # z0, background cAMP signal

x0=0.01*np.ones(N)
y0=0.08*np.ones(N)
z0=0
x_trace=np.zeros((N,nSteps)); x_trace[:,0] = x0
y_trace=np.zeros((N,nSteps)); y_trace[:,0] = y0
z_trace=np.zeros((nSteps,1)); z_trace[0] = z0
Kamino_pop = Kamino2017_pop_SCnoise(x0,y0,z0, KaminoPopParam)

for i in range(nSteps-1):
    x_next,y_next,z_next=Kamino_pop.update(cAMPext_influx_trace[i],dt)
    x_trace[:,i+1] = x_next
    y_trace[:,i+1] = y_next
    z_trace[i+1] = z_next
         
later_portion = 0.25 # start count peaks after this X total simulation time
y_trace_norm =  (y_trace-Nh_Kamino_offset)/Nh_Kamino # height normalization
y_trace_norm_later=y_trace_norm[:,math.floor(nSteps * later_portion):] # the later part of trace
y_trace_norm_mean = np.mean(y_trace_norm,axis=0)
y_trace_norm_mean_later=y_trace_norm_mean[math.floor(nSteps * later_portion):] # the later part of trace

t_plot_Kamino = np.array(t)/Nt_Kamino
t_plot_Kamino_short = t_plot_Kamino[0:math.floor(nSteps * (1-later_portion))] 

pop_max = np.amax(y_trace_norm_mean_later); pop_min = np.amin(y_trace_norm_mean_later); 
pk_find_thresh = 0.02
PkPos, PkProperties = find_peaks(y_trace_norm_mean_later, prominence=pk_find_thresh)

#pop_max = np.amax(y_trace_norm_mean_later); pop_min = np.amin(y_trace_norm_mean_later); 
#pk_find_prm = 0.4
#PkPos, PkProperties = find_peaks(y_trace_norm_mean_later, prominence=pk_find_prm*(pop_max-pop_min))

if len(PkPos) == 0:
    firing_rate = 0; height = 0
else: 
    firing_rate = len(PkPos)/(t_tot*(1-later_portion))
    height = np.mean(PkProperties["prominences"])

# Check find_peaks
title = 'Kamino 2017, noise=' +str(sigma)+', density='+str(rho)+', J='+str(gamma)+'\n ,pk_find_thresh='+str(pk_find_thresh)
fig = plt.figure()
plt.plot(y_trace_norm_mean_later)
plt.plot(PkPos, y_trace_norm_mean_later[PkPos], "x")
plt.title(title+' FR= '+str(firing_rate))

# plot population mean and selected single cell traces
SC_traces_idx = [0,2,4,6,8,10]
# title = r'Martiel 1987, $cAMP_{e}$ input = '+str(cAMPext_influx)
plot_POP_oscillation(t_plot_Kamino_short,y_trace_norm_mean_later,cAMPext_influx,
                     t_tot*(1-later_portion), 0,title, y_trace_norm_later, SC_traces_idx)

##  check simulation traces
#fig,ax = plt.subplots()
#t_plot_Kamino = np.array(t)/Nt_Kamino
#liney = ax.plot(t_plot_Kamino,y_trace, label= r'y, $cAMP_{cyto}$')
#linez = ax.plot(t_plot_Kamino,z_trace, label =  r'z,$cAMP_{ext}$')
#ax.set_ylim([-0.05,0.4])
#ax.set_xlabel('Time')
#ax.set_ylabel('x,y,z')
#ax.set_title('Kamino 2017 group oscillation, rho= '+str(rho))
#leg = ax.legend()
#ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
#plt.show()

# #%% Get the oscillation period 
# y_trace=np.array(y_trace) # convert list to array
# later_portion = 0.2 # start count peaks after this X total simulation time
# y_trace_later=y_trace[math.floor(len(t)*later_portion):] # the later part of trace
# PkPos, PkProperties = find_peaks(y_trace_later, prominence=(0.02,100))

# ## Check find_peaks
# #plt.plot(y_trace_later) 
# #plt.plot(PkPos, y_trace_later[PkPos], "x")

# Kamino_pop_osc_period = (1-later_portion)*t_tot / len(PkPos)
# print('group oscillation period for Kamino is '+str(Kamino_pop_osc_period))

#%% Sgro 2015
from Sgro2015_agent_and_pop_FUN import Sgro2015_pop

from Params import SgroPopParam
sigma = 0.15 # noise strength
N = 100 # number of cells in the population
rho = 10**(-3.5); j = 0.5

A0=-1.5*np.ones(N) # + np.random.uniform(-sigma,sigma,N); ###########
R0=-0.5*np.ones(N) # + np.random.uniform(-sigma,sigma,N); 
cAMPext0 = 0
Sgro_pop=Sgro2015_pop(A0,R0,cAMPext0, SgroPopParam)

dt=0.005 ; t_tot=20; t=list(np.arange(0,t_tot*Nt_Sgro,dt))
nSteps = len(t)
# define extracellular stim trace
alphafval =0
time_separation = 0

# initializations
A_trace_orig=np.zeros((N,len(t)))
R_trace_orig=np.zeros((N,len(t)))
cAMPext_trace = np.zeros((len(t),1))
A_trace_orig[:,0]=A0 
R_trace_orig[:,0]=R0
cAMPext_trace[0] = cAMPext0

for i in range(len(t)-1):
    A_now=A_trace_orig[:,i]
    R_now=R_trace_orig[:,i]
    cAMPext_now = cAMPext_trace[i]
    
    A_next,R_next, cAMPext_next = Sgro_pop.update( dt, time_separation,alphafval)
    A_trace_orig[:,i] = A_next
    R_trace_orig[:,i] = R_next
    cAMPext_trace[i] = cAMPext_next

    
# Traces
later_portion = 0.25
A_trace_norm=(A_trace_orig - Nh_Sgro_offset)/Nh_Sgro;
A_trace_norm_later = A_trace_norm[:,math.floor(nSteps * later_portion):] # the later part of trace
A_trace_norm_mean = np.mean(A_trace_norm,axis=0)
A_trace_norm_mean_later=A_trace_norm_mean[math.floor(nSteps * later_portion):] 
# cAMPext_trace_plot = np.array(cAMPext_trace)

t_plot_Sgro = np.array(t)/Nt_Sgro
t_plot_Sgro_short = t_plot_Sgro[0:math.floor(nSteps * (1-later_portion))] 

pk_find_thresh = 0.4
PkPos, PkProperties = find_peaks(A_trace_norm_mean_later, prominence=pk_find_thresh)

if len(PkPos) == 0:
    firing_rate = 0; height = 0
else: 
    firing_rate = len(PkPos)/(t_tot*(1-later_portion))
    height = np.mean(PkProperties["prominences"])

# Check find_peaks
title = 'Sgro2015, noise=' +str(sigma)+', density='+str(rho)+', J='+str(j)+'\n ,pk_find_thresh='+str(pk_find_thresh)
fig = plt.figure()
plt.plot(A_trace_norm_mean_later)
plt.plot(PkPos, A_trace_norm_mean_later[PkPos], "x")
plt.title(title+' FR= '+str(firing_rate))

# plot population mean and selected single cell traces
SC_traces_idx = [0,2,4,6,8,10]
# title = r'Martiel 1987, $cAMP_{e}$ input = '+str(cAMPext_influx)
plot_POP_oscillation(t_plot_Sgro_short,A_trace_norm_mean_later,0,
                     t_tot*(1-later_portion), 0,title, A_trace_norm_later, SC_traces_idx)


# #  check simulation traces
# label_font_size=25; trace_width=3; tick_font_size=18

# fig,ax = plt.subplots()
# ax.plot(t_plot_Sgro,A_trace_mean_plot,linewidth=trace_width)
# # ax.plot(t_plot_Sgro,A_trace_plot[2,:],linewidth=trace_width, label= r'activator, $cAMP_{cyto}$')
# # ax.set_ylim([-0.2,1.3])
# ax.set_xlabel('Time')
# ax.set_ylabel(r'Activator, $cAMP_{cyto}$')
# ax.set_title(r'Sgro 2015 group oscillation, $rho$= '+str(rho)+', j= '+str(j)+', time separation= '+str(time_separation))
# ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
# plt.show()

## Get the extracellular cAMP oscillation height
#later_portion = 0.2 # start count peaks after this X total simulation time
#cAMPext_trace_plot_later=cAMPext_trace_plot[math.floor(len(t)*later_portion):] # the later part of trace
#PkPos, PkProperties = find_peaks(cAMPext_trace_plot_later[:,0], prominence=(0.5,500))
#plt.figure()
#plt.plot(cAMPext_trace_plot_later)
#plt.plot(PkPos, cAMPext_trace_plot_later[PkPos], "x")

## Get the oscillation period
#later_portion = 0.2 # start count peaks after this X total simulation time
#A_trace_mean_plot_later=A_trace_mean_plot[math.floor(len(t)*later_portion):] # the later part of trace
#PkPos, PkProperties = find_peaks(A_trace_mean_plot_later, prominence=(0.5,1.5))
## Check find_peaks
#plt.plot(A_trace_mean_plot_later)
#plt.plot(PkPos, A_trace_mean_plot_later[PkPos], "x")
#
#Sgro_pop_osc_period = (1-later_portion)*t_tot / len(PkPos)
#print('group oscillation period for Sgro is '+str(Sgro_pop_osc_period))

#%% Gregor 2010

from Gregor2010_agent_and_pop_FUN import Gregor2010_pop

from Params import GregorPopParam
Amax=20;  Abas=0.4 # uM
Nc=100 # Num of cells
eta=0.02 # noise stength

ext_input = 0
time_separation = 0


campCyto0 = 7.5*np.ones(Nc)
sinthetai0 = (campCyto0*2-Amax-Abas)/(-Amax+Abas)
thetai0 = np.arcsin(sinthetai0)
campExt0 = 0 # Vc*St/Sc*rho/K*c_sec*1/Nc*np.sum(campCyto0);

Gregor_pop=Gregor2010_pop(campCyto0, thetai0, campExt0, GregorPopParam)

dt=0.005; t_tot=20; t=list(np.arange(0,t_tot*Nt_Gregor,dt))
nSteps = len(t)
# initializations
gregor_thetai_trace=np.zeros((Nc,len(t))) 
gregor_campCyto_trace=np.zeros((Nc,len(t))) 
gregor_campExt_trace=np.zeros(len(t)) 

for i in range(nSteps - 1):
    thetai_now=gregor_thetai_trace[:,i]
    campCyto_now=gregor_campCyto_trace[:,i]
    campExt_now=gregor_campExt_trace[i]
    thetai_next, campCyto_next, campExt_next = Gregor_pop.update(dt,eta,rho,k,Vt,time_separation,ext_input)
    gregor_thetai_trace[:,i+1] = thetai_next
    gregor_campCyto_trace[:,i+1] = campCyto_next
    gregor_campExt_trace[i+1] = campExt_next
    
#Traces
# gregor_thetai_trace= np.array(gregor_thetai_trace) 
# gregor_campExt_trace = np.array(gregor_campExt_trace)
later_portion = 0.25
gregor_campCyto_trace_norm= np.array(gregor_campCyto_trace)/Nh_Gregor
gregor_campCyto_trace_norm_later=gregor_campCyto_trace_norm[:,math.floor(nSteps * later_portion):] # the later part of trace
gregor_campCyto_trace_norm_mean = np.mean(gregor_campCyto_trace_norm,axis=0)
gregor_campCyto_trace_norm_mean_later=gregor_campCyto_trace_norm_mean[math.floor(nSteps * later_portion):]

t_plot_Gregor = np.array(t)/Nt_Gregor
t_plot_Gregor_short = t_plot_Gregor[0:math.floor(nSteps * (1-later_portion))] 

pop_max = np.amax(y_trace_norm_mean_later); pop_min = np.amin(y_trace_norm_mean_later); 
pk_find_thresh = 0.5
PkPos, PkProperties = find_peaks(gregor_campCyto_trace_norm_mean_later, prominence=pk_find_thresh)

#pop_max = np.amax(y_trace_norm_mean_later); pop_min = np.amin(y_trace_norm_mean_later); 
#pk_find_prm = 0.4
#PkPos, PkProperties = find_peaks(y_trace_norm_mean_later, prominence=pk_find_prm*(pop_max-pop_min))

if len(PkPos) == 0:
    firing_rate = 0; height = 0
else: 
    firing_rate = len(PkPos)/(t_tot*(1-later_portion))
    height = np.mean(PkProperties["prominences"])

# Check find_peaks
title = 'Gregor 2010, noise=' +str(eta)+', density='+str(rho)+', J='+str(k)+'\n ,pk_find_thresh='+str(pk_find_thresh)
fig = plt.figure()
plt.plot(gregor_campCyto_trace_norm_mean_later)
plt.plot(PkPos, gregor_campCyto_trace_norm_mean_later[PkPos], "x")
plt.title(title+' FR= '+str(firing_rate))

# plot population mean and selected single cell traces
SC_traces_idx = [0,2,4,6,8,10]
# title = r'Martiel 1987, $cAMP_{e}$ input = '+str(cAMPext_influx)
plot_POP_oscillation(t_plot_Gregor_short,gregor_campCyto_trace_norm_mean_later,cAMPext_influx,
                     t_tot*(1-later_portion), 0,title, gregor_campCyto_trace_norm_later, SC_traces_idx)


# #  check simulation traces
# label_font_size=25; trace_width=3; tick_font_size=18

# fig,ax = plt.subplots()
# ax.plot(t_plot_Gregor,gregor_campCyto_trace_mean,linewidth=trace_width, label= r'activator, $cAMP_{cyto}$')
# # ax.plot(t_plot_Gregor,gregor_campExt_trace,linewidth=trace_width, label= r'$cAMP_{ext}$')
# # ax.plot(t_plot_Sgro,A_trace_plot[2,:],linewidth=trace_width, label= r'activator, $cAMP_{cyto}$')
# # ax.set_ylim([-0.2,1.3])
# ax.set_xlabel('Time')
# ax.set_ylabel('Activator')
# ax.set_title(r'Gregor 2010 group oscillation, $rho$= '+str(rho)+', k= '+str(k)+', time separation= '+str(time_separation))
# leg = ax.legend()
# ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
# plt.show()

##Get the oscillation period
#later_portion = 0.2 # start count peaks after this X total simulation time
#gregor_campCyto_trace_mean_later=gregor_campCyto_trace_mean[math.floor(len(t)*later_portion):] # the later part of trace
#PkPos, PkProperties = find_peaks(gregor_campCyto_trace_mean_later, prominence=(2,25))

## Check find_peaks
#plt.plot(gregor_campCyto_trace_mean_later)
#plt.plot(PkPos, gregor_campCyto_trace_mean_later[PkPos], "x")

#Gregor_pop_osc_period = (1-later_portion)*t_tot / len(PkPos)
#print('group oscillation period for Gregor is '+str(Gregor_pop_osc_period))
#%%
label_font_size=30
trace_width=5
tick_font_size=20
mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']
fig5 = plt.figure(figsize=(12, 12))
grid = plt.GridSpec(3, 1, wspace=0.2, hspace=0.2)

ax2= fig5.add_subplot(grid[:, 0])
ax2.plot(t_plot_Goldbeter, b_trace, color=mycolors[0],linewidth=trace_width, label='Goldbeter 1987')
ax2.plot(t_plot_Maeda_short,cAMPi_trace_later,color=mycolors[1], linewidth=trace_width, label='Maeda & Loomis 1998')
ax2.plot(t_plot_Gregor,gregor_campCyto_trace_mean,color=mycolors[2],linewidth=trace_width, label='Gregor 2010')
ax2.plot(t_plot_Sgro, A_trace_mean_plot,color=mycolors[5],linewidth=trace_width, label='Sgro & Mehta 2015')
ax2.plot(t_plot_Kamino, y_trace, color=mycolors[7],linewidth=trace_width, label='Kamino & Sawai 2017')

ax2.set_xlim([-1,16])
ax2.set_ylim([-0.2,1.3])
ax2.set_ylabel(r'$cAMP_{cyto}$, A.U.',fontsize=label_font_size)
ax2.set_xlabel('Time, A.U.',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
leg = ax2.legend(frameon=False,loc='upper center',ncol=2,prop={'size': 17});
plt.show()

#%% Save all outputs in npz file
np.savez('pop_oscillation_SCnoise_200409.npz', 
         t_plot_Goldbeter_short = t_plot_Goldbeter_short , 
         b_trace_norm_later = b_trace_norm_later, b_trace_norm_mean_later = b_trace_norm_mean_later,
         t_plot_Maeda_short=t_plot_Maeda_short, 
         cAMPi_trace_norm_later = cAMPi_trace_norm_later,
         cAMPi_trace_norm_mean_later = cAMPi_trace_norm_mean_later,
         t_plot_Gregor_short=t_plot_Gregor_short, 
         gregor_campCyto_trace_norm_later =gregor_campCyto_trace_norm_later,
         gregor_campCyto_trace_norm_mean_later = gregor_campCyto_trace_norm_mean_later,
         t_plot_Sgro_short=t_plot_Sgro_short, 
         A_trace_norm_later = A_trace_norm_later,
         A_trace_norm_mean_later = A_trace_norm_mean_later,
         t_plot_Kamino_short=t_plot_Kamino_short, 
         y_trace_norm_later = y_trace_norm_later,
         y_trace_norm_mean_later = y_trace_norm_mean_later)
         
#%% experimental data
my_dir = r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/'
Sgro2015Figure6excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure6')

#%% load saved npz output file
npzfile = np.load('pop_oscillation_SCnoise_200409.npz')

t_plot_Goldbeter_short =  npzfile['t_plot_Goldbeter_short'] ; 
b_trace_norm_later = npzfile['b_trace_norm_later']
b_trace_norm_mean_later = npzfile['b_trace_norm_mean_later']

t_plot_Maeda_short=npzfile['t_plot_Maeda_short'] ; 
cAMPi_trace_norm_later=npzfile['cAMPi_trace_norm_later']
cAMPi_trace_norm_mean_later = npzfile['cAMPi_trace_norm_mean_later']

t_plot_Gregor_short=npzfile['t_plot_Gregor_short']; 
gregor_campCyto_trace_norm_later = npzfile['gregor_campCyto_trace_norm_later']
gregor_campCyto_trace_norm_mean_later = npzfile['gregor_campCyto_trace_norm_mean_later']

t_plot_Sgro_short=npzfile['t_plot_Sgro_short']; 
A_trace_norm_later = npzfile['A_trace_norm_later']
A_trace_norm_mean_later = npzfile['A_trace_norm_mean_later']

t_plot_Kamino_short= npzfile['t_plot_Kamino_short']; 
y_trace_norm_later = npzfile['y_trace_norm_later']
y_trace_norm_mean_later = npzfile['y_trace_norm_mean_later']
#%% Plot all models in separate subplots
#title_font_size = 22
#label_font_size = 22
#tick_font_size = 16
#legend_font_size = 12
#trace_width = 3

abcd_font_size = 28
label_font_size=24
title_font_size = 26
sublabel_font_size = 22
trace_width=3
tick_font_size=20
SC_traces_idx = [0,2,4,6,8,10]

mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']
#fig3 = plt.figure(figsize=(16, 14))
#grid = plt.GridSpec(2, 2, wspace=0.15, hspace=0.25)
fig3 = plt.figure(figsize=(12,10))
grid = plt.GridSpec(3, 2, wspace=0.2, hspace=0.8)

ax0= fig3.add_subplot(grid[0, 0])
ax0.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["No External cAMP Mean Trace"],
                              color = 'k',linewidth=trace_width)
ax0.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax0.set_title('Experiment', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax0.set_ylim([-0.1,0.5]); ax0.set_xlim([0,120])
ax0.set_xlabel('Time (min)', size=sublabel_font_size)
ax0.text(-0.08 , 1.2, 'A', ha='center',va='center',
     transform = ax0.transAxes, color = 'b', fontsize=abcd_font_size)

ax1= fig3.add_subplot(grid[0, 1])
for this_idx in SC_traces_idx:
    this_trace = b_trace_norm_later[this_idx,:] 
    ax1.plot(t_plot_Goldbeter_short,this_trace, color='k',alpha=0.6, linewidth=2)     
# Plot population mean
ax1.plot(t_plot_Goldbeter_short,b_trace_norm_mean_later,
         color=mycolors[0],alpha=0.8, linewidth=trace_width)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Martiel 1987', color=mycolors[0], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax1.set_ylim([-0.3,1.7]); ax1.set_xlim([0,15])
ax1.text(-0.08 , 1.2, 'B', ha='center',va='center',
     transform = ax1.transAxes, color = 'g', fontsize=abcd_font_size)

ax2= fig3.add_subplot(grid[1,0])
for this_idx in SC_traces_idx:
    this_trace = cAMPi_trace_norm_later[this_idx,:] 
    ax2.plot(t_plot_Maeda_short,this_trace, color='k',alpha=0.6, linewidth=2)     
# Plot population mean
ax2.plot(t_plot_Maeda_short,cAMPi_trace_norm_mean_later,
         color=mycolors[1],alpha=0.8, linewidth=trace_width)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Maeda & Loomis 2004',color=mycolors[1], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax2.set_ylim([0.1,0.6]); ax2.set_xlim([0,15])
ax2.text(-0.08 , 1.2, 'C', ha='center',va='center',
     transform = ax2.transAxes, color = 'g', fontsize=abcd_font_size)

ax3= fig3.add_subplot(grid[1, 1])
for this_idx in SC_traces_idx:
    this_trace = gregor_campCyto_trace_norm_later[this_idx,:] 
    ax3.plot(t_plot_Gregor_short,this_trace, color='k',alpha=0.6, linewidth=2)     
# Plot population mean
ax3.plot(t_plot_Gregor_short,gregor_campCyto_trace_norm_mean_later,
         color=mycolors[2],alpha=0.8, linewidth=trace_width)
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Gregor 2010', color=mycolors[2], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax3.set_ylim([-0.2,1.5]); ax3.set_xlim([0,15])
ax3.text(-0.08 , 1.2, 'D', ha='center',va='center',
     transform = ax3.transAxes, color = 'g', fontsize=abcd_font_size)

ax4= fig3.add_subplot(grid[2, 0])
for this_idx in SC_traces_idx:
    this_trace = A_trace_norm_later[this_idx,:] 
    ax4.plot(t_plot_Sgro_short,this_trace, color='k',alpha=0.6, linewidth=2)     
# Plot population mean
ax4.plot(t_plot_Sgro_short,A_trace_norm_mean_later,
         color=mycolors[5],alpha=0.8, linewidth=trace_width)
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('Sgro 2015', color=mycolors[5], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax4.set_ylim([-0.3,1.3]); ax4.set_xlim([0,15])
ax4.text(-0.08 , 1.2, 'E', ha='center',va='center',
     transform = ax4.transAxes, color = 'g', fontsize=abcd_font_size)

ax5= fig3.add_subplot(grid[2,1])
for this_idx in SC_traces_idx:
    this_trace = y_trace_norm_later[this_idx,:] 
    ax5.plot(t_plot_Kamino_short,this_trace, color='k',alpha=0.6, linewidth=2)     
# Plot population mean
ax5.plot(t_plot_Kamino_short,y_trace_norm_mean_later,
         color=mycolors[7],alpha=0.8, linewidth=trace_width)
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5.set_title('Kamino 2017', color=mycolors[7], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax5.set_ylim([-0.5,1.5]); ax5.set_xlim([0,15])
ax5.text(-0.08 , 1.2, 'F', ha='center',va='center',
     transform = ax5.transAxes, color = 'g', fontsize=abcd_font_size)

fig3.text(0.04, 0.5, r'$cAMP_{i}$', fontsize=label_font_size,va='center', rotation='vertical')
fig3.text(0.5, 0.04, 'Time, A.U.', fontsize=label_font_size, ha='center')
plt.show()
