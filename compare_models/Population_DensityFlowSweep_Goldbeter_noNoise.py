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

# import code and data locations
from config import CODE_DIR, output_dir
import Goldbeter1987Func as gbf

# import normalization parameters
with open(CODE_DIR / 'normParams.pickle','rb') as f:
    Nt_Gregor,Nt_Sgro,Nt_Goldbeter,Nt_Maeda,Nt_Kamino,Nh_Gregor,Nh_Sgro,Nh_Goldbeter,Nh_Maeda,Nh_Kamino,Nh_Kamino_offset,Nh_Sgro_offset = pickle.load(f)

# define variable arrays 
kcVals = np.linspace(0,100,25)
hVals = 1/np.logspace(-2,1,25)

# get index values from job submission 
kcIndex = int(sys.argv[1])
hIndex = int(sys.argv[2])

# set model parameters
k1 = 0.036      # per min
k2 = 0.666      # per min
L1 = 10
L2 = 0.005 
c = 10          # 0.15 ~ 50
lamda =0.01
theta =0.01
e =  1          # 0.108 # compared to 1
q=4000
sig = 0.6       # 0.57 # compared to 0.6
v = 12
k= 4            # k prime in the paper
ki = 1.7        # 0.958 # compared to 1.7 
kt = 0.9
kc = kcVals[kcIndex]        # 3.58 # compared to 5.4
h = hVals[hIndex]           # ratio of intracellular to extracellular volume, rho 
N = 100         # number of cells in the population
Goldbeter3PopParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h, 'N':N}

# set initial values 
p0=0.8
a0=3
b0=0.9
g0=0

# set simulation time parameters 
dt=0.001
t_tot=30*Nt_Goldbeter
t=np.arange(0,t_tot,dt)

# define extracellular stim amounts 
cAMPe_trace = np.zeros(len(t))

# initialize random generator and set random values 
rng = np.random.default_rng(seed=1)
sigma = 0      # noise strength
# the run draws its noise from rng at each time step: math.sqrt(dt)*rng.normal(0,1,N)*sigma

# initialize population
pop = gbf.Population_scNoise([p0,a0,b0,g0],Goldbeter3PopParam,t)

# run population 
pop.run(dt,rng,a0,cAMPe_trace,sigma=sigma)

# normalize back 80% of trace 
traceStart = 0.2
bs = pop.b_mean[int(len(t)*traceStart):]/Nh_Goldbeter

# analyze trace for firing rate and peak height 
firingRate = 0
peakHeight = 0
PkPos, PkProperties = find_peaks(bs, prominence = 5/Nh_Goldbeter)
if len(PkPos) != 0:
    firingRate = len(PkPos)/(t_tot/Nt_Goldbeter*(1-traceStart))
    peakHeight = np.mean(PkProperties["prominences"])

ts = t[int(len(t)*traceStart):]/Nt_Goldbeter

# export to pickle file
with open(output_dir('pop_firingRate_data/no_noise/Goldbeter') /
          ('kc_%.2i_h_%.2i.pickle' % (kcIndex,hIndex)),'wb') as f:
    pickle.dump([ts,bs,firingRate,peakHeight],f,pickle.HIGHEST_PROTOCOL)

