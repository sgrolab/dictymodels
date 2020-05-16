# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 2019

@author: Chuqiao

Population collective oscillations

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

# Normalization parameters
from Params import NormParams
for key,val in NormParams.items():
        exec(key + '=val')

#%% Sgro 2015
from Sgro2015_agent_and_pop_FUN import Sgro2015_pop

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5; sigma = 0.15 # noise strength
N = 100 # number of cells in the population
rho = 10**(-3.5); j = 0.5
SgroPopParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':sigma,'N':N,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'flux_thrs':0, 'rho': rho,'j': j}

A0=-1.5*np.ones(N) # + np.random.uniform(-sigma,sigma,N); ###########
R0=-0.5*np.ones(N) # + np.random.uniform(-sigma,sigma,N); 
cAMPext0 = 0
Sgro_pop=Sgro2015_pop(A0,R0,cAMPext0, SgroPopParam)

dt=0.005 ; t_tot=15*Nt_Sgro; t=list(np.arange(0,t_tot,dt))

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
A_trace_offset=1.5
A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
A_trace_plot=(A_trace_orig - Nh_Sgro_offset)/Nh_Sgro;
A_trace_mean_plot = np.mean(A_trace_plot,axis = 0)
t_plot_Sgro = np.array(t)/Nt_Sgro

cAMPext_trace_plot = np.array(cAMPext_trace)

# Convert into np array

#  check simulation traces
label_font_size=25; trace_width=3; tick_font_size=18

fig,ax = plt.subplots()
ax.plot(t_plot_Sgro,A_trace_mean_plot,linewidth=trace_width)
# ax.plot(t_plot_Sgro,A_trace_plot[2,:],linewidth=trace_width, label= r'activator, $cAMP_{cyto}$')
# ax.set_ylim([-0.2,1.3])
ax.set_xlabel('Time')
ax.set_ylabel(r'Activator, $cAMP_{cyto}$')
ax.set_title(r'Sgro 2015 group oscillation, $rho$= '+str(rho)+', j= '+str(j)+', time separation= '+str(time_separation))
ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
plt.show()

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

Amax=20;  Abas=0.4 # uM
w=2*pi/6 # min-1
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
ext_input = 0
time_separation = 0

GregorPopParam={'Amax':Amax,'Abas':Abas,'w':w,'Vc':Vc,'St':St,'Sc':Sc,'K':K,\
            'c_sec':c_sec,'c_excite':c_excite,'Nc':Nc}

campCyto0 = 7.5*np.ones(Nc)
sinthetai0 = (campCyto0*2-Amax-Abas)/(-Amax+Abas)
thetai0 = np.arcsin(sinthetai0)
campExt0 = 0 # Vc*St/Sc*rho/K*c_sec*1/Nc*np.sum(campCyto0);

Gregor_pop=Gregor2010_pop(campCyto0, thetai0, campExt0, GregorPopParam)

dt=0.005; t_tot=15*Nt_Gregor; t=list(np.arange(0,t_tot,dt))

eta=0.002 # noise stength
# initializations
gregor_thetai_trace=np.zeros((Nc,len(t))) 
gregor_campCyto_trace=np.zeros((Nc,len(t))) 
gregor_campExt_trace=np.zeros(len(t)) 

for i in range(len(t)-1):
    thetai_now=gregor_thetai_trace[:,i]
    campCyto_now=gregor_campCyto_trace[:,i]
    campExt_now=gregor_campExt_trace[i]
    thetai_next, campCyto_next, campExt_next = Gregor_pop.update(dt,eta,rho,k,Vt,time_separation,ext_input)
    gregor_thetai_trace[:,i+1] = thetai_next
    gregor_campCyto_trace[:,i+1] = campCyto_next
    gregor_campExt_trace[i+1] = campExt_next
    
#Traces
gregor_thetai_trace= np.array(gregor_thetai_trace) 
gregor_campCyto_trace= np.array(gregor_campCyto_trace)/Nh_Gregor
gregor_campCyto_trace_mean= np.mean(gregor_campCyto_trace,axis = 0)
gregor_campExt_trace = np.array(gregor_campExt_trace)
t_plot_Gregor = np.array(t)/Nt_Gregor


#  check simulation traces
label_font_size=25; trace_width=3; tick_font_size=18

fig,ax = plt.subplots()
ax.plot(t_plot_Gregor,gregor_campCyto_trace_mean,linewidth=trace_width, label= r'activator, $cAMP_{cyto}$')
# ax.plot(t_plot_Gregor,gregor_campExt_trace,linewidth=trace_width, label= r'$cAMP_{ext}$')
# ax.plot(t_plot_Sgro,A_trace_plot[2,:],linewidth=trace_width, label= r'activator, $cAMP_{cyto}$')
# ax.set_ylim([-0.2,1.3])
ax.set_xlabel('Time')
ax.set_ylabel('Activator')
ax.set_title(r'Gregor 2010 group oscillation, $rho$= '+str(rho)+', k= '+str(k)+', time separation= '+str(time_separation))
leg = ax.legend()
ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
plt.show()

##Get the oscillation period
#later_portion = 0.2 # start count peaks after this X total simulation time
#gregor_campCyto_trace_mean_later=gregor_campCyto_trace_mean[math.floor(len(t)*later_portion):] # the later part of trace
#PkPos, PkProperties = find_peaks(gregor_campCyto_trace_mean_later, prominence=(2,25))

## Check find_peaks
#plt.plot(gregor_campCyto_trace_mean_later)
#plt.plot(PkPos, gregor_campCyto_trace_mean_later[PkPos], "x")

#Gregor_pop_osc_period = (1-later_portion)*t_tot / len(PkPos)
#print('group oscillation period for Gregor is '+str(Gregor_pop_osc_period))

#%% Golbeter 1987,  Table II/ Fig 3 parameters
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

Goldbeter3PopParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}

p0=0.8; a0=3; b0=0.9; g0=0

# Fig 2, 4 variable model, autonomous oscillation
Goldbeter3_pop = Goldbeter1987_pop_3var([1,1],[p0,a0,b0,g0],Goldbeter3PopParam)
# initializations
p_trace=[p0]; b_trace=[b0]; g_trace=[g0]

dt=0.001; t_tot=15*Nt_Goldbeter; t=list(np.arange(0,t_tot,dt))
campExt_influx = 0
for i in range(len(t)-1):
    p_next,b_next,g_next= Goldbeter3_pop.update(dt,a0,campExt_influx)
    p_trace.append(p_next)
    b_trace.append(b_next)
    g_trace.append(g_next)
        
   
# Convert into np array
b_trace = np.array(b_trace);  b_trace = b_trace/Nh_Goldbeter
p_trace = np.array(p_trace);  p_trace = p_trace/np.amax(p_trace)
t_plot_Goldbeter = np.array(t)/Nt_Goldbeter

#  check simulation traces
label_font_size=25; trace_width=3; tick_font_size=18

fig,ax = plt.subplots()
ax.plot(t_plot_Goldbeter,b_trace,linewidth=trace_width, label= r'b, $cAMP_{cyto}$')
# ax.plot(t_plot_Goldbeter,p_trace, linewidth=trace_width,label = r'p, $R_{act}/R_{tot}$')
# ax.set_ylim([-0.2,1.3])
ax.set_xlabel('Time')
ax.set_ylabel('b, p')
ax.set_title('Goldber 1987 group oscillation, fig 3 parameters')
leg = ax.legend()
ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
plt.show()

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

#%% Laub & Loomis 1998
#from LaubLoomis1998_agent_and_pop_FUN import MaedaLoomis2003_agent
## parameters from Maeda & Loomis 2003 paper
#k1=2.0; k2=0.9; k3=2.5; k4=1.5; k5=0.6
#k6=0.8; k7=1.0; k8=1.3; k9=0.3; k10=0.8
#k11=0.7; k12=4.9; k13=23; k14=4.5
#MaedaAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
#            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
#            'k13':k13,'k14':k14}
### parameters from Maeda & Loomis 2003 paper
##k1=1.4; k2=0.9; k3=2.5; k4=1.5; k5=0.6
##k6=0.8; k7=2.0; k8=1.3; k9=0.7; k10=1.0
##k11=0.3; k12=3.1; k13=1.8; k14=1.5
##LaubAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
##            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
##            'k13':k13,'k14':k14}
#ACA0=0.1; PKA0=0.1; ERK20=0.2; RegA0=0.1; cAMPi0=0.5; 
#cAMPe0=0.05; CAR10=0.1
#state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]
#MaedaLoomis_pop=MaedaLoomis2003_agent([1,1],state0,MaedaAgentParam)
#
## initializations
#ACA_trace=[ACA0]; PKA_trace=[PKA0]; ERK2_trace=[ERK20]
#RegA_trace=[RegA0]; cAMPi_trace=[cAMPi0]; # cAMPe_trace=[cAMPe0]
#CAR1_trace=[CAR10]
#
#dt=0.0005; t_tot=500; t=list(np.arange(0,t_tot,dt))
#signal_input = 'none'
#
#for i in range(len(t)-1):
#    ACA_now=ACA_trace[i]
#    PKA_now=PKA_trace[i]
#    ERK2_now=ERK2_trace[i]
#    RegA_now=RegA_trace[i]
#    cAMPi_now=cAMPi_trace[i]
#    cAMPe_now=cAMPi_trace[i]
#    CAR1_now=CAR1_trace[i]
#    
#    ACA_next,PKA_next,ERK2_next,RegA_next,\
#    cAMPi_next,cAMPe_next,CAR1_next=MaedaLoomis_pop.update(dt,signal_input)
#    
#    ACA_trace.append(ACA_next)
#    PKA_trace.append(PKA_next)
#    ERK2_trace.append(ERK2_next)
#    RegA_trace.append(RegA_next)
#    cAMPi_trace.append(cAMPi_next)
#    # cAMPe_trace.append(cAMPe_next)
#    CAR1_trace.append(CAR1_next)
#    
#
#ERK2_trace = np.array(ERK2_trace) ; ERK2_trace = ERK2_trace/np.amax(ERK2_trace)
#cAMPi_trace = np.array(cAMPi_trace); cAMPi_trace = cAMPi_trace/np.amax(cAMPi_trace)
#t_plot_Maeda = np.array(t)
#
##  check simulation traces
#label_font_size=25; trace_width=3; tick_font_size=18
#
#fig,ax = plt.subplots()
#t_plot_Kamino = np.array(t)
#ax.plot(t_plot_Maeda,cAMPi_trace,linewidth=trace_width, label= r'$cAMP_{cyto}$')
#ax.plot(t_plot_Maeda,ERK2_trace, linewidth=trace_width,label = r'ERK2')
## ax.set_ylim([-0.2,1.3])
#ax.set_xlabel('Time')
#ax.set_ylabel(r'$cAMP_{cyto}$, ERK2')
#ax.set_title('Maeda Loomis 2003 group oscillation')
#leg = ax.legend()
#ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
#plt.show()
#
## Get the oscillation period 
#cAMPi_trace=np.array(cAMPi_trace) # convert list to array
#later_portion = 0.2 # start count peaks after this X total simulation time
#cAMPi_trace_later=cAMPi_trace[math.floor(len(t)*later_portion):] # the later part of trace
#PkPos, PkProperties = find_peaks(cAMPi_trace_later, prominence=(0,1.2))
#
### Check find_peaks
##plt.plot(cAMPi_trace_later)
##plt.plot(PkPos, cAMPi_trace_later[PkPos], "x")
#
#Maeda_pop_osc_period = (1-later_portion)*t_tot / len(PkPos)
#print('group oscillation period for Maeda & Loomis is '+str(Maeda_pop_osc_period))

#%% Maeda & Loomis 1998
from MaedaLoomis2004_agent_and_pop_FUN import MaedaLoomis2004_pop
# parameters from Maeda & Loomis 2004 paper
k1=2.0; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=1.0; k8=1.3; k9=0.3; k10=0.8
k11=0.7; k12=4.9; k13=23; k14=4.5
MaedaAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}
## parameters from Maeda & Loomis 2003 paper
#k1=1.4; k2=0.9; k3=2.5; k4=1.5; k5=0.6
#k6=0.8; k7=2.0; k8=1.3; k9=0.7; k10=1.0
#k11=0.3; k12=3.1; k13=1.8; k14=1.5
#LaubAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
#            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
#            'k13':k13,'k14':k14}
ACA0=0.1; PKA0=0.1; ERK20=0.1; RegA0=0.1; cAMPi0=0.01; 
cAMPe0=0.1; CAR10=0.1
state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]
MaedaLoomis_pop=MaedaLoomis2004_pop([1,1],state0,MaedaAgentParam)

# initializations
ACA_trace=[ACA0]; PKA_trace=[PKA0]; ERK2_trace=[ERK20]
RegA_trace=[RegA0]; cAMPi_trace=[cAMPi0]; # cAMPe_trace=[cAMPe0]
CAR1_trace=[CAR10]

dt=0.001; t_tot=60*Nt_Maeda; t=list(np.arange(0,t_tot,dt))
signal_input = 0

for i in range(len(t)-1):
    ACA_now=ACA_trace[i]
    PKA_now=PKA_trace[i]
    ERK2_now=ERK2_trace[i]
    RegA_now=RegA_trace[i]
    cAMPi_now=cAMPi_trace[i]
    cAMPe_now=cAMPi_trace[i]
    CAR1_now=CAR1_trace[i]
    
    ACA_next,PKA_next,ERK2_next,RegA_next,\
    cAMPi_next,cAMPe_next,CAR1_next=MaedaLoomis_pop.update(dt,signal_input)
    
    ACA_trace.append(ACA_next)
    PKA_trace.append(PKA_next)
    ERK2_trace.append(ERK2_next)
    RegA_trace.append(RegA_next)
    cAMPi_trace.append(cAMPi_next)
    # cAMPe_trace.append(cAMPe_next)
    CAR1_trace.append(CAR1_next)
    

ERK2_trace = np.array(ERK2_trace) ; ERK2_trace = ERK2_trace/np.amax(ERK2_trace)
cAMPi_trace = np.array(cAMPi_trace)/Nh_Maeda
t_plot_Maeda = np.array(t)/Nt_Maeda
t_plot_Maeda_short = t_plot_Maeda[0:int(len(t)*0.25)]
cAMPi_trace_later = cAMPi_trace[0:int(len(t)*0.25)] #cAMPi_trace[int(len(t_plot_Maeda)*0.75):]

#  check simulation traces
label_font_size=25; trace_width=3; tick_font_size=18

fig,ax = plt.subplots()
ax.plot(t_plot_Maeda_short,cAMPi_trace_later,linewidth=trace_width, label= r'$cAMP_{cyto}$')
# ax.plot(t_plot_Maeda,ERK2_trace, linewidth=trace_width,label = r'ERK2')
# ax.set_ylim([-0.2,1.3])
ax.set_xlabel('Time')
ax.set_ylabel(r'$cAMP_{cyto}$, ERK2')
ax.set_title('Maeda Loomis 2003 group oscillation')
leg = ax.legend()
ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
plt.show()

## Get the oscillation period 
#cAMPi_trace=np.array(cAMPi_trace) # convert list to array
#later_portion = 0.2 # start count peaks after this X total simulation time
#cAMPi_trace_later=cAMPi_trace[math.floor(len(t)*later_portion):] # the later part of trace
#PkPos, PkProperties = find_peaks(cAMPi_trace_later, prominence=(0,1.2))
#
### Check find_peaks
##plt.plot(cAMPi_trace_later)
##plt.plot(PkPos, cAMPi_trace_later[PkPos], "x")
#
#Maeda_pop_osc_period = (1-later_portion)*t_tot / len(PkPos)
#print('group oscillation period for Maeda & Loomis is '+str(Maeda_pop_osc_period))
#%% Kamino 2017, fig 5D group oscillations
from Kamino2017_agent_and_pop_FUN import Kamino2017_pop
   
tau=1.5; n=2; K=4; kt=2; delta=0
gamma = 3; rho = 1 # cell density


dt=0.001
t_tot= 15*Nt_Kamino                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
t=np.arange(0,t_tot,dt)
z0_influx_trace=np.zeros(len(t)) # z0, background cAMP signal

# start simulation
tic = perf_counter() 

Param={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
       'gamma':gamma,'rho':rho}

x0=0.01
y0=0.08
z0=0.01
x_trace=[x0]
y_trace=[y0]
z_trace=[z0]

test_agent=Kamino2017_pop([x0,y0,z0],Param)

for i in range(len(t)-1):
    x_now=x_trace[i]
    y_now=y_trace[i]
    z0_influx = 0
    
    x_next,y_next,z_next=test_agent.update(z0_influx,dt)
    x_trace.append(x_next)
    y_trace.append(y_next)
    z_trace.append(z_next)
toc = perf_counter()  # show how long have passed
print('time passed:'+str(toc-tic))    
y_trace = np.array(y_trace)/Nh_Kamino
#  check simulation traces
fig,ax = plt.subplots()
t_plot_Kamino = np.array(t)/Nt_Kamino
liney = ax.plot(t_plot_Kamino,y_trace, label= r'y, $cAMP_{cyto}$')
linez = ax.plot(t_plot_Kamino,z_trace, label =  r'z,$cAMP_{ext}$')
# ax.set_ylim([-0.05,0.4])
ax.set_xlabel('Time')
ax.set_ylabel('x,y,z')
ax.set_title('Kamino 2017 group oscillation, rho= '+str(rho))
leg = ax.legend()
ax.legend( frameon=False,loc='upper center',ncol=2,prop={'size': 15})
plt.show()

#%% Get the oscillation period 
y_trace=np.array(y_trace) # convert list to array
later_portion = 0.2 # start count peaks after this X total simulation time
y_trace_later=y_trace[math.floor(len(t)*later_portion):] # the later part of trace
PkPos, PkProperties = find_peaks(y_trace_later, prominence=(0.02,100))

## Check find_peaks
#plt.plot(y_trace_later)
#plt.plot(PkPos, y_trace_later[PkPos], "x")

Kamino_pop_osc_period = (1-later_portion)*t_tot / len(PkPos)
print('group oscillation period for Kamino is '+str(Kamino_pop_osc_period))

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
np.savez('pop_oscillation_200225.npz', 
         t_plot_Goldbeter = t_plot_Goldbeter , b_trace=b_trace,
         t_plot_Maeda_short=t_plot_Maeda_short, cAMPi_trace_later=cAMPi_trace_later,
         t_plot_Gregor=t_plot_Gregor, gregor_campCyto_trace_mean=gregor_campCyto_trace_mean,
         t_plot_Sgro=t_plot_Sgro, A_trace_mean_plot=A_trace_mean_plot,
         t_plot_Kamino=t_plot_Kamino, y_trace=y_trace)
         
#%% experimental data
my_dir = r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/'
Sgro2015Figure6excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure6')

#%% load saved npz output file
npzfile = np.load('pop_oscillation_200225.npz')
t_plot_Goldbeter =  npzfile['t_plot_Goldbeter'] ; b_trace = npzfile['b_trace']
t_plot_Maeda_short=npzfile['t_plot_Maeda_short'] ; cAMPi_trace_later=npzfile['cAMPi_trace_later']
t_plot_Gregor=npzfile['t_plot_Gregor']; gregor_campCyto_trace_mean=npzfile['gregor_campCyto_trace_mean']
t_plot_Sgro=npzfile['t_plot_Sgro']; A_trace_mean_plot=npzfile['A_trace_mean_plot']
t_plot_Kamino= npzfile['t_plot_Kamino']; y_trace=npzfile['y_trace']
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
ax1.plot(t_plot_Goldbeter,b_trace,color=mycolors[0],linewidth=trace_width)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Martiel 1987', color=mycolors[0], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax1.set_ylim([-0.2,1.6]); ax1.set_xlim([0,15])
ax1.text(-0.08 , 1.2, 'B', ha='center',va='center',
     transform = ax1.transAxes, color = 'g', fontsize=abcd_font_size)


ax2= fig3.add_subplot(grid[1,0])
ax2.plot(t_plot_Maeda_short,cAMPi_trace_later,color=mycolors[1],linewidth=trace_width)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Maeda & Loomis 2004',color=mycolors[1], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax2.set_ylim([0.1,0.6]); ax2.set_xlim([0,15])
ax2.text(-0.08 , 1.2, 'C', ha='center',va='center',
     transform = ax2.transAxes, color = 'g', fontsize=abcd_font_size)

ax3= fig3.add_subplot(grid[1, 1])
ax3.plot(t_plot_Gregor,gregor_campCyto_trace_mean,color=mycolors[2],linewidth=trace_width)
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Gregor 2010', color=mycolors[2], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax3.set_ylim([-0.2,1.2]); ax3.set_xlim([0,15])
ax3.text(-0.08 , 1.2, 'D', ha='center',va='center',
     transform = ax3.transAxes, color = 'g', fontsize=abcd_font_size)

ax4= fig3.add_subplot(grid[2, 0])
ax4.plot(t_plot_Sgro,A_trace_mean_plot,color=mycolors[5],linewidth=trace_width)
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('Sgro 2015', color=mycolors[5], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax4.set_ylim([-0.2,1.2]); ax4.set_xlim([0,15])
ax4.text(-0.08 , 1.2, 'E', ha='center',va='center',
     transform = ax4.transAxes, color = 'g', fontsize=abcd_font_size)

ax5= fig3.add_subplot(grid[2,1])
ax5.plot(t_plot_Kamino, y_trace, color=mycolors[7],linewidth=trace_width)
ax5.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax5.set_title('Kamino 2017', color=mycolors[7], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax5.set_ylim([-0.2,1.2]); ax5.set_xlim([0,15])
ax5.text(-0.08 , 1.2, 'F', ha='center',va='center',
     transform = ax5.transAxes, color = 'g', fontsize=abcd_font_size)

fig3.text(0.04, 0.5, r'$cAMP_{i}$', fontsize=label_font_size,va='center', rotation='vertical')
fig3.text(0.5, 0.04, 'Time, A.U.', fontsize=label_font_size, ha='center')
plt.show()
