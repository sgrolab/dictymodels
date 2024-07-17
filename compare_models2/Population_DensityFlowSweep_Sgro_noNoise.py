# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 2019

@author: Chuqiao

Population oscillations with added [cAMP]ext

"""
import sys
import pickle
import numpy as np
import math

from scipy.signal import find_peaks

# add path for function files
sys.path.insert(0,'//groups/sgro/sgrolab/mark/dicty_proj/dictymodels/compare_models2')
import Sgro2015Func as sg

# import normalization parameters
with open('//groups/sgro/sgrolab/mark/dicty_proj/dictymodels/compare_models2/normParams.pickle','rb') as f:
    Nt_Gregor,Nt_Sgro,Nt_Goldbeter,Nt_Maeda,Nt_Kamino,Nh_Gregor,Nh_Sgro,Nh_Goldbeter,Nh_Maeda,Nh_Kamino,Nh_Kamino_offset,Nh_Sgro_offset = pickle.load(f)

# define variable arrays 
jVals = np.linspace(0,1,26) # 21 j_arr=np.array([0.5]) # 
rhoVals = np.logspace(-6,-3,26) # 26 rho_arr=np.array([10**(-3.5)]) # 

# get index values from job submission 
jIndex = int(sys.argv[1])
rhoIndex = int(sys.argv[2])

# set model parameters 
e=0.1
tauA=0.09
tauR=tauA/e
g=0.5
sigma = 0        # noise strength
rho = rhoVals[rhoIndex]
j = jVals[jIndex]
N = 100             # number of cells in the population
SgroPopParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g, 'c0':1.2, \
              'sigma':sigma,'N':N,'a':0.058,'alpha0':800,'alpha_pde':1000,\
              'Kd':1e-5,'S':1e6,'flux_thrs':0, 'rho':rho, 'j':j}

# set initial values
A0 = -1.5
R0 = -0.5
cAMPext0 = 0
    
# set simulation time parameters 
dt=0.005
t_tot=25*Nt_Sgro
t = np.arange(0,t_tot,dt)

# set time separation variable 
time_separation = 0

# initialize random generator and set random values 
rng = np.random.default_rng(seed=1)
r = math.sqrt(dt)*rng.normal(0,1,size = (len(t),N))

# define extracellular stim amounts
extracellularA_trace = np.zeros(len(t))

# initialize cell population 
Sgro_pop = sg.Sgro2015_pop(A0,R0,cAMPext0,SgroPopParam,t)

# run population 
Sgro_pop.run(dt,time_separation,extracellularA_trace,r)

# normalize back 80% of trace 
traceStart = 0.2
As = (np.mean(Sgro_pop.A,axis=0)[int(len(t)*traceStart):]-Nh_Sgro_offset)/Nh_Sgro

# analyze trace for firing rate and peak height 
firingRate = 0
peakHeight = 0
PkPos, PkProperties = find_peaks(As, prominence=(0.4,2))
if len(PkPos) != 0:
    firingRate = len(PkPos)/(t_tot/Nt_Sgro*(1-traceStart))
    peakHeight = np.mean(PkProperties["prominences"])

# save adjusted time trace
ts = t[int(len(t)*traceStart):]/Nt_Sgro

with open('//groups/sgro/sgrolab/mark/dicty_proj/pop_firingRate_data/no_noise/Sgro/j_%.2i_rho_%.2i.pickle' % (jIndex,rhoIndex),'wb') as f:
    pickle.dump([ts,As,firingRate,peakHeight],f,pickle.HIGHEST_PROTOCOL)




