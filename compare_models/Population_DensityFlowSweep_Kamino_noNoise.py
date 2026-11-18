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
from popHistory import find_peaks    # scipy.signal's is slow on a noisy trace: see there

# import data locations
from config import output_dir
import Kamino2017Func as kf

# import normalization parameters
from normParams import Nt_Kamino, Nh_Kamino, Nh_Kamino_offset

# define variable arrays 
gammaVals = np.linspace(0,100,31)
rhoVals = np.logspace(0,2,31)

# get index values from job submission 
gammaIndex = int(sys.argv[1])
rhoIndex = int(sys.argv[2])

# set model parameters
tau=1.5
n=2
K=4
kt=2
delta = 0.01
gamma = gammaVals[gammaIndex]   # dilution rate 
rho = rhoVals[rhoIndex]         # cell density
N = 100     # number of cells in the population 
params={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
       'N':N, 'gamma':gamma,'rho':rho}

# set initial values
x0=0.01
y0=0.08
z0=0

# set simulation time parameters
dt=0.0001
t_tot = 30 * Nt_Kamino
t = np.arange(0,t_tot,dt)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              

# initialize random generator and set random values 
rng = np.random.default_rng(seed=1)
sigma =  0   # noise strength
# the run draws its noise from rng at each time step: math.sqrt(dt)*rng.normal(0,1,N)*sigma

# define cAMPe amounts  
camp_input_trace=np.zeros(len(t))

# initialize cell population 
pop = kf.Population_scNoise([x0,y0,z0],params,t)

# run simulation
pop.run(camp_input_trace,dt,rng,sigma=sigma)

# normalize back 80% of trace 
traceStart = 0.2
As = (pop.y_mean[int(len(t)*traceStart):]-Nh_Kamino_offset)/Nh_Kamino

# analyze trace for firing rate and peak height 
firingRate = 0
peakHeight = 0
PkPos, PkProperties = find_peaks(As, prominence=(0.02,100))
if len(PkPos) != 0:
    firingRate = len(PkPos)/(t_tot/Nt_Kamino*(1-traceStart))
    peakHeight = np.mean(PkProperties["prominences"])

# save adjusted time trace
ts = t[int(len(t)*traceStart):]/Nt_Kamino

with open(output_dir('pop_firingRate_data/no_noise/Kamino') /
          ('j_%.2i_rho_%.2i.pickle' % (gammaIndex,rhoIndex)),'wb') as f:
    pickle.dump([ts,As,firingRate,peakHeight],f,pickle.HIGHEST_PROTOCOL)
