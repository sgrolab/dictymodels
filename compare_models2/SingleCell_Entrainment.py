#%% Single Cell, Step vs Ramp Input

import os, sys, pickle 
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 

# add path for function files
sys.path.insert(0,'//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

# import normalization parameters
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/normParams.pickle','rb') as f:
    Nt_Gregor,Nt_Sgro,Nt_Goldbeter,Nt_Maeda,Nt_Kamino,Nh_Gregor,Nh_Sgro,Nh_Goldbeter,Nh_Maeda,Nh_Kamino,Nh_Kamino_offset,Nh_Sgro_offset = pickle.load(f)

#%% Sgro 

import Sgro2015Func as sf

# set parameters 
SgroAgentParam={'e':0.1,'tauA':0.09,'tauR':0.9,'g':0.5,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}

# set time parameters
dt=0.001

# set initial values 
A0=-1.5
R0=-0.5

# set noise amount 
rng = np.random.default_rng(seed=1)

# set number of cycles 
nCycles = 8

# set cAMPe input wave period and width and amplitidue  
peakPeriods = np.linspace(0.8,1.8,8)
peakWidths = np.linspace(0.3,1.4,10)
cAMPe = 1

# preallocate entrainment array 
Rs = np.zeros([len(peakWidths),len(peakPeriods)])

for i in range(len(peakWidths)):
    for j in range(len(peakPeriods)):
        if peakWidths[i] >= peakPeriods[j]:
            Rs[i,j] = np.nan
        else:
            
            # define cAMPe 
            stim_1cycle = np.zeros(int(peakPeriods[j]*Nt_Sgro/dt))
            stim_1cycle[0:int(peakWidths[i]*Nt_Sgro/dt)] = cAMPe 
            signal_trace = np.concatenate((np.zeros(int(1*Nt_Sgro/dt)), np.tile(stim_1cycle,(nCycles)),np.zeros(int(1*Nt_Sgro/dt))),axis=0);
            
            # define noise 
            r = np.sqrt(dt) * rng.normal(0,1,size = (len(signal_trace)))

            # initialize cell 
            cell = sf.Cell([1,1],[A0,R0],SgroAgentParam,signal_trace)
            
            # run simulation 
            cell.run(dt,signal_trace,r)

            # normalize cAMPi trace 
            A = (cell.A-Nh_Sgro_offset)/Nh_Sgro
            
            # get initial peak 
            initialPeak = A[int((1+0*peakPeriods[j])*Nt_Sgro/dt):int((1+1*peakPeriods[j])*Nt_Sgro/dt)]
            
            # preallocate correlation coefficient array 
            sampleRs = np.zeros(nCycles-1)
            
            # find correlation between first peak and each subsequent one 
            for k in range(1,nCycles):
                
                # calculate time point at which peak starts 
                t_peakStart = int((1+k*peakPeriods[j])*Nt_Sgro/dt)
                
                followPeak = A[t_peakStart:t_peakStart+len(initialPeak)]
                correlation = np.corrcoef(initialPeak,followPeak)
                sampleRs[k-1] = correlation[1,0]
            
            # store mean of correlation coefficients 
            Rs[i,j] = np.mean(sampleRs)

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_Entrainment_data/entrainment_Sgro.pickle','wb') as f:
    pickle.dump(Rs,f,pickle.HIGHEST_PROTOCOL)



#%% Goldbeter 

import Goldbeter1987Func as gbf

# set parameters 
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

# set time parameters
dt=0.001
    
# set initial values 
p0=0.8
a0=3
b0=0.9
g0=0

# set number of cycles 
nCycles = 8

# set cAMPe input wave period and width and amplitidue  
peakPeriods = np.linspace(0.8,1.8,8)
peakWidths = np.linspace(0.3,1.4,10)
cAMPe = 1

# preallocate entrainment array 
Rs = np.zeros([len(peakWidths),len(peakPeriods)])

for i in range(len(peakWidths)):
    for j in range(len(peakPeriods)):
        if peakWidths[i] >= peakPeriods[j]:
            Rs[i,j] = np.nan
        else:
            
            # define cAMPe 
            stim_1cycle = np.zeros(int(peakPeriods[j]*Nt_Goldbeter/dt))
            stim_1cycle[0:int(peakWidths[i]*Nt_Goldbeter/dt)] = cAMPe 
            signal_trace = np.concatenate((np.zeros(int(1*Nt_Goldbeter/dt)), np.tile(stim_1cycle,(nCycles)),np.zeros(int(1*Nt_Goldbeter/dt))),axis=0);
            
            # initialize cell 
            cell = gbf.Cell(Goldbeter3AgentParam,[p0,b0,g0],signal_trace)
            
            # run cell 
            cell.run(dt,a0,signal_trace)
            
            # normalize cAMPi trace 
            b = cell.b/Nh_Goldbeter
            
            # get initial peak 
            initialPeak = b[int((1+0*peakPeriods[j])*Nt_Goldbeter/dt):int((1+1*peakPeriods[j])*Nt_Goldbeter/dt)]
            
            # preallocate correlation coefficient array 
            sampleRs = np.zeros(nCycles-1)
            
            # find correlation between first peak and each subsequent one 
            for k in range(1,nCycles):
                
                # calculate time point at which peak starts 
                t_peakStart = int((1+k*peakPeriods[j])*Nt_Goldbeter/dt)
                
                followPeak = b[t_peakStart:t_peakStart+len(initialPeak)]
                correlation = np.corrcoef(initialPeak,followPeak)
                sampleRs[k-1] = correlation[1,0]
            
            # store mean of correlation coefficients 
            Rs[i,j] = np.mean(sampleRs)

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_Entrainment_data/entrainment_Goldbeter.pickle','wb') as f:
    pickle.dump(Rs,f,pickle.HIGHEST_PROTOCOL)


#%% Maeda Loomis

import MaedaLoomis2004Func as mlf

# set parameters
k1=2.0; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=1.0; k8=1.3; k9=0.3; k10=0.8
k11=0.7; k12=4.9; k13=23; k14=4.5
MaedaAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}

# set time parameters
dt=0.0005

# set initial values 
ACA0=0.1
PKA0=0.1
ERK20=0.1
RegA0=0.1
cAMPi0=0.01
cAMPe0=0.1
CAR10=0.1

# set number of cycles 
nCycles = 8

# set cAMPe input wave period and width and amplitidue  
peakPeriods = np.linspace(0.8,1.8,8)
peakWidths = np.linspace(0.3,1.4,10)
cAMPe = 1

# preallocate entrainment array 
Rs = np.zeros([len(peakWidths),len(peakPeriods)])

for i in range(len(peakWidths)):
    for j in range(len(peakPeriods)):
        if peakWidths[i] >= peakPeriods[j]:
            Rs[i,j] = np.nan
        else:
            
            # define cAMPe 
            stim_1cycle = np.zeros(int(peakPeriods[j]*Nt_Maeda/dt))
            stim_1cycle[0:int(peakWidths[i]*Nt_Maeda/dt)] = cAMPe 
            signal_trace = np.concatenate((np.zeros(int(1*Nt_Maeda/dt)), np.tile(stim_1cycle,(nCycles)),np.zeros(int(1*Nt_Maeda/dt))),axis=0);
            
            # initialize cell 
            cell = mlf.Cell([0,0],[ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10],MaedaAgentParam,signal_trace)
            
            # run cell
            cell.run(dt,signal_trace)
            
            # normalize cAMPi trace 
            cAMPi = cell.cAMPi/Nh_Maeda
            
            # get initial peak 
            initialPeak = cAMPi[int((1+0*peakPeriods[j])*Nt_Maeda/dt):int((1+1*peakPeriods[j])*Nt_Maeda/dt)]
            
            # preallocate correlation coefficient array 
            sampleRs = np.zeros(nCycles-1)
            
            # find correlation between first peak and each subsequent one 
            for k in range(1,nCycles):
                
                # calculate time point at which peak starts 
                t_peakStart = int((1+k*peakPeriods[j])*Nt_Maeda/dt)
                
                followPeak = cAMPi[t_peakStart:t_peakStart+len(initialPeak)]
                correlation = np.corrcoef(initialPeak,followPeak)
                sampleRs[k-1] = correlation[1,0]
            
            # store mean of correlation coefficients 
            Rs[i,j] = np.mean(sampleRs)

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_Entrainment_data/entrainment_Maeda.pickle','wb') as f:
    pickle.dump(Rs,f,pickle.HIGHEST_PROTOCOL)


#%% Kamino 

import Kamino2017Func as kf

# set parameters 
tau=1.5
n=2
K=4
kt=2
delta=0.01
gamma=3
rho= 0.01
KaminoAgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}

# set time parameters 
dt=0.001

# set initial values 
x0=0.01
y0=0.06
z0=0.005

# set number of cycles 
nCycles = 8

# set cAMPe input wave period and width and amplitidue  
peakPeriods = np.linspace(0.8,1.8,8)
peakWidths = np.linspace(0.3,1.4,10)
cAMPe = 1

# preallocate entrainment array 
Rs = np.zeros([len(peakWidths),len(peakPeriods)])

for i in range(len(peakWidths)):
    for j in range(len(peakPeriods)):
        if peakWidths[i] >= peakPeriods[j]:
            Rs[i,j] = np.nan
        else:
            
            # define cAMPe 
            stim_1cycle = np.zeros(int(peakPeriods[j]*Nt_Kamino/dt))
            stim_1cycle[0:int(peakWidths[i]*Nt_Kamino/dt)] = cAMPe 
            signal_trace = np.concatenate((np.zeros(int(1*Nt_Kamino/dt)), np.tile(stim_1cycle,(nCycles)),np.zeros(int(1*Nt_Kamino/dt))),axis=0);
            
            # initialize cell 
            cell = kf.Cell([x0,y0,z0],KaminoAgentParam,signal_trace)
            
            # run cell
            cell.run(dt,signal_trace)

            # normalize cAMPi trace 
            y = (cell.y-Nh_Kamino_offset)/Nh_Kamino
            
            # get initial peak 
            initialPeak = y[int((1+0*peakPeriods[j])*Nt_Kamino/dt):int((1+1*peakPeriods[j])*Nt_Kamino/dt)]
            
            # preallocate correlation coefficient array 
            sampleRs = np.zeros(nCycles-1)
            
            # find correlation between first peak and each subsequent one 
            for k in range(1,nCycles):
                
                # calculate time point at which peak starts 
                t_peakStart = int((1+k*peakPeriods[j])*Nt_Kamino/dt)
                
                followPeak = y[t_peakStart:t_peakStart+len(initialPeak)]
                correlation = np.corrcoef(initialPeak,followPeak)
                sampleRs[k-1] = correlation[1,0]
            
            # store mean of correlation coefficients 
            Rs[i,j] = np.mean(sampleRs)

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_Entrainment_data/entrainment_Kamino.pickle','wb') as f:
    pickle.dump(Rs,f,pickle.HIGHEST_PROTOCOL)

