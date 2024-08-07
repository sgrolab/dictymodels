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
import MaedaLoomis2004Func as mlf

# import normalization parameters
with open('//groups/sgro/sgrolab/mark/dicty_proj/dictymodels/compare_models2/normParams.pickle','rb') as f:
    Nt_Gregor,Nt_Sgro,Nt_Goldbeter,Nt_Maeda,Nt_Kamino,Nh_Gregor,Nh_Sgro,Nh_Goldbeter,Nh_Maeda,Nh_Kamino,Nh_Kamino_offset,Nh_Sgro_offset = pickle.load(f)

# define variable arrays 
gammaVals = np.linspace(0,100,26)
rhoVals = np.logspace(-0.1,2.7,26)

# get index values from job submission 
gammaIndex = int(sys.argv[1])
rhoIndex = int(sys.argv[2])

# parameters from Maeda & Loomis 2004 paper
k1=2.0; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=1.0; k8=1.3; k9=0.3; k10=0.8
k11=0.7; k12=4.9; k13=23; k14=4.5
N = 100 # number of cells in a population

MaedaPopParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14, 'N':N}

# set initial values 
ACA0=0.1
PKA0=0.1
ERK20=0.1
RegA0=0.1
cAMPi0=0.01 
cAMPe0=0.1
CAR10=0.1
initialVals = [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]

# set simulation time parameters
dt = 0.0001
t_tot = 60*Nt_Maeda
t = np.arange(0,t_tot,dt)

# initialize random generator and set random values 
rng = np.random.default_rng(seed=1)
sigma = 0      # noise strength
r = math.sqrt(dt)*rng.normal(0,1,size = (len(t),N)) * sigma

# define extracellular cAMP input amounts 
cAMPe_in = np.zeros(len(t))

# initialize population 
pop = mlf.Population_scNoise([1,1],initialVals,MaedaPopParam,t)

# run simulation 
pop.run(dt,cAMPe_in,r,rhoVals[rhoIndex],gammaVals[gammaIndex])

# normalize back 50% of cAMPi trace 
traceStart = 0.5
cAMPi = np.mean(pop.cAMPi,axis=0)[int(len(t)*traceStart):]/Nh_Maeda

# analyze trace for firing rate and peak height 
firingRate = 0
peakHeight = 0
pk_find_prm = 0.7
PkPos, PkProperties = find_peaks(cAMPi,prominence=((np.amax(cAMPi)-np.amin(cAMPi))*pk_find_prm,np.amax(cAMPi)))
if len(PkPos) != 0:
    firingRate = len(PkPos)/(t_tot/Nt_Maeda*(1-traceStart))
    peakHeight = np.mean(PkProperties["prominences"])

# normalize back 50% of time values 
ts = t[int(len(t)*traceStart):]/Nt_Maeda

# export time and cAMPi data
with open('//groups/sgro/sgrolab/mark/dicty_proj/pop_firingRate_data/no_noise/Maeda/g_%.2i_rho_%.2i.pickle' % (gammaIndex,rhoIndex),'wb') as f:
    pickle.dump([ts,cAMPi,firingRate,peakHeight],f,pickle.HIGHEST_PROTOCOL)
