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
import Gregor2010Func as gf

# import normalization parameters
with open('//groups/sgro/sgrolab/mark/dicty_proj/dictymodels/compare_models2/normParams.pickle','rb') as f:
    Nt_Gregor,Nt_Sgro,Nt_Goldbeter,Nt_Maeda,Nt_Kamino,Nh_Gregor,Nh_Sgro,Nh_Goldbeter,Nh_Maeda,Nh_Kamino,Nh_Kamino_offset,Nh_Sgro_offset = pickle.load(f)

# define variable arrays 
rhoVals = np.logspace(-3.5,1,26) 
kVals = np.linspace(1,50,26) # np.array([25]) #

# get index values from job submission 
rhoIndex = int(sys.argv[1])
kIndex = int(sys.argv[2])

# set model parameters
Amax=20                 # uM
Abas=0.4                # uM
w=2*math.pi/6           # min-1
Vc=1.1e-9               # ml
St=1.33                 # cm2
Sc=1.3e-6               # cm2
K=0.0004                # uM, 400 pM
c_sec= 3.6              # min-1
c_excite=1.01           # min-1
Nc=100                  # Num of cells
rho = rhoVals[rhoIndex] # 1/ml
Vt = 1                  # chamber size (ml)
k = kVals[kIndex]       # ml/min
GregorPopParam={'Amax':Amax,'Abas':Abas,'w':w,'Vc':Vc,'St':St,'Sc':Sc,'K':K,\
            'c_sec':c_sec,'c_excite':c_excite,'Nc':Nc, 'rho':rho, 'Vt':Vt,'k':k}
 
# set initial values
campCyto0 = 7.5*np.ones(Nc)
thetai0 = np.arcsin((campCyto0*2-Amax-Abas)/(-Amax+Abas))
campExt0 = 0 

# set simulation time parameters
dt = 0.001
t_tot = 25*Nt_Gregor
t = np.arange(0,t_tot,dt)

# define extracellular stim amounts 
extracellularA_trace = np.zeros(len(t))

# set time separation variable 
time_separation = 0

# initialize random generator and set random steps 
eta = 0.002 # noise stength
rng = np.random.default_rng(seed=1)
r = math.sqrt(dt)*rng.normal(0,1,size = (len(t),Nc))

# initialize cell population 
Gregor_pop = gf.Gregor2010_pop(campCyto0, thetai0, campExt0, GregorPopParam, t)

# run population
Gregor_pop.run(dt,eta,rho,k,Vt,time_separation,extracellularA_trace,r)

# normalize back 80% of cAMPi trace 
traceStart = 0.2
cAMPi = np.mean(Gregor_pop.cAMPi,axis=0)[int(len(t)*traceStart):]/Nh_Gregor

# analyze trace for firing rate and peak height 
firingRate = 0
peakHeight = 0
PkPos, PkProperties = find_peaks(cAMPi,prominence=(0.5,30))
if len(PkPos) != 0:
    firingRate = len(PkPos)/(t_tot/Nt_Gregor*(1-traceStart))
    peakHeight = np.mean(PkProperties["prominences"])

# normalize back 50% of time values 
ts = t[int(len(t)*traceStart):]/Nt_Gregor

# export time and cAMPi data
with open('//groups/sgro/sgrolab/mark/dicty_proj/pop_firingRate_data/noise/Gregor/rho_%.2i_k_%.2i.pickle' % (rhoIndex,kIndex),'wb') as f:
    pickle.dump([ts,cAMPi,firingRate,peakHeight],f,pickle.HIGHEST_PROTOCOL)
