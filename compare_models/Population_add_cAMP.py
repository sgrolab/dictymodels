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

# add path for function files
sys.path.insert(0,'//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

# import normalization parameters
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/normParams.pickle','rb') as f:
    Nt_Gregor,Nt_Sgro,Nt_Goldbeter,Nt_Maeda,Nt_Kamino,Nh_Gregor,Nh_Sgro,Nh_Goldbeter,Nh_Maeda,Nh_Kamino,Nh_Kamino_offset,Nh_Sgro_offset = pickle.load(f)

#%% Sgro 2015

import Sgro2015Func as sg

# set model parameters 
e=0.1
tauA=0.09
tauR=tauA/e
g=0.5
sigma = 0.15        # noise strength
rho = 10**(-3.5)
j = 0.5
N = 100             # number of cells in the population
SgroPopParam={'e':0.1,'tauA':0.09,'tauR':tauR,'g':0.5, 'c0':1.2, \
              'sigma':0.15,'N':N,'a':0.058,'alpha0':800,'alpha_pde':1000,\
              'Kd':1e-5,'S':1e6,'flux_thrs':0, 'rho':rho, 'j':j}

# set initial values
A0 = -1.5
R0 = -0.5
cAMPext0 = 0
    
# set simulation time parameters 
dt=0.005
t_tot=30*Nt_Sgro
t = np.arange(0,t_tot,dt)

# define extracellular stim amounts
extracellularAs = np.array([10,20,100])

# set time separation variable 
time_separation = 0

# initialize random generator and set random values 
rng = np.random.default_rng(seed=1)
r = math.sqrt(dt)*rng.normal(0,1,size = (len(t),N))

# preallocate A trace arrays
As = np.zeros((len(extracellularAs),len(t)))
As_sc = np.zeros((len(extracellularAs),N,len(t)))

# iterate through each extracellular A amount 
for i in range(len(extracellularAs)):
    
    # define extracellular A trace with input applied halfway
    extracellularA_trace = np.zeros(len(t))
    extracellularA_trace[int(len(t)/2):] = extracellularAs[i]
    
    # initialize cell population 
    Sgro_pop = sg.Sgro2015_pop(A0,R0,cAMPext0,SgroPopParam,t)
    
    # run population 
    Sgro_pop.run(dt,time_separation,extracellularA_trace,r)
    
    # Save A traces to arrays
    A_trace_offset=1.5
    As_sc[i,:,:] = (Sgro_pop.A+A_trace_offset)/Nh_Sgro
    As[i,:] = np.mean(As_sc[i,:,:], axis = 0)

# save adjusted time trace
ts = t/Nt_Sgro

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/popAddcAMP_data/popAddcAMP_Sgro.pickle','wb') as f:
    pickle.dump([ts,As,As_sc],f,pickle.HIGHEST_PROTOCOL)


#%% Gregor 2010

import Gregor2010Func as gf

# set model parameters
Amax=20;        # uM
Abas=0.4        # uM
w=2*math.pi/6   # min-1
Vc=1.1e-9       # ml
St=1.33         # cm2
Sc=1.3e-6       # cm2
K=0.0004        # uM, 400 pM
c_sec= 3.6      # min-1
c_excite=1.01   # min-1
Nc=100           # Num of cells
rho = 1/12      # 1/ml
Vt = 1          # chamber size (ml)
k = 5           # ml/min
GregorPopParam={'Amax':Amax,'Abas':Abas,'w':w,'Vc':Vc,'St':St,'Sc':Sc,'K':K,\
            'c_sec':c_sec,'c_excite':c_excite,'Nc':Nc, 'rho':rho, 'Vt':Vt,'k':k}

# set initial values
campCyto0 = 7.5*np.ones(Nc)
thetai0 = np.arcsin((campCyto0*2-Amax-Abas)/(-Amax+Abas))
campExt0 = 0 # Vc*St/Sc*rho/K*c_sec*1/Nc*np.sum(campCyto0);

# set simulation time parameters
dt = 0.001
t_tot = 30*Nt_Gregor
t = np.arange(0,t_tot,dt)
eta = 0.002 # noise stength

# define extracellular stim amounts 
ext_input_arr =  np.array([0.001,10,1000])

# set time separation variable 
time_separation = 0

# initialize random generator and set random steps 
rng = np.random.default_rng(seed=1)
r = math.sqrt(dt)*rng.normal(0,1,size = (len(t),N))

# fix random seed
np.random.seed(1)
r = math.sqrt(dt)*np.random.normal(0,1,size = (len(t),Nc))

# preallocate trace arrays
thetais = np.zeros([len(ext_input_arr),Nc,len(t)])
cAMPes = np.zeros([len(ext_input_arr),len(t)])
cAMPis = np.zeros((len(ext_input_arr),len(t)))
cAMPis_sc = np.zeros((len(ext_input_arr),Nc,len(t)))

# iterate through each extracellular cAMP amount 
for i in range(len(ext_input_arr)):

    # define extracellular A trace with input applied halfway
    extracellularA_trace = np.zeros(len(t))
    extracellularA_trace[int(len(t)/2):] = ext_input_arr[i]
    
    # initialize cell population 
    Gregor_pop = gf.Gregor2010_pop(campCyto0, thetai0, campExt0, GregorPopParam, t)
    
    # run population
    Gregor_pop.run(dt,eta,rho,k,Vt,time_separation,extracellularA_trace,r)
    
    # save traces 
    thetais[i] = Gregor_pop.thetai
    cAMPes[i] = Gregor_pop.cAMPe
    cAMPis_sc[i] = Gregor_pop.cAMPi/Nh_Gregor   # adjust trace
    cAMPis[i] = np.mean(cAMPis_sc[i],axis=0)

# save adjusted time trace 
ts = np.array(t)/Nt_Gregor

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/popAddcAMP_data/popAddcAMP_Gregor.pickle','wb') as f:
    pickle.dump([ts,thetais,cAMPes,cAMPis_sc,cAMPis],f,pickle.HIGHEST_PROTOCOL)


#%% Golbeter 1987, Table II/ Fig 3 parameters

import Goldbeter1987Func as gbf

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
kc = 5.4        # 3.58 # compared to 5.4
h = 5           # ratio of intracellular to extracellular volume, rho 
Goldbeter3AgentParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}

# set initial values 
p0=0.8
a0=3
b0=0.9
g0=0

# set simulation time parameters 
dt=0.0005
t_tot=30*Nt_Goldbeter
t=np.arange(0,t_tot,dt)

# define extracellular stim amounts 
cAMPe_amounts = np.array([0.01,0.02,0.1])

# preallocate arrays
ps = np.zeros([len(cAMPe_amounts),len(t)])
bs = np.zeros_like(ps)
gs = np.zeros_like(ps)

# run model for each cAMPe amount 
for i in range(len(cAMPe_amounts)):
    
    # define cAMPe trace 
    cAMPe_trace = np.zeros(len(t))
    cAMPe_trace[int(len(t)/2):] = cAMPe_amounts[i]
    
    # initialize population
    pop = gbf.Population([p0,a0,b0,g0],Goldbeter3AgentParam,t)
    
    # run population 
    pop.run(dt,a0,cAMPe_trace)
    
    # save traces
    ps[i] = pop.p
    bs[i] = pop.b
    gs[i] = pop.g

bs = bs/Nh_Goldbeter
ts = t//Nt_Goldbeter

# export to pickle file
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/popAddcAMP_data/popAddcAMP_Goldbeter.pickle','wb') as f:
    pickle.dump([ts,ps,bs,gs],f,pickle.HIGHEST_PROTOCOL)

#%% Maeda & Loomis 2004
import MaedaLoomis2004Func as mlf

# parameters from Maeda & Loomis 2004 paper
k1=2.0; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=1.0; k8=1.3; k9=0.3; k10=0.8
k11=0.7; k12=4.9; k13=23; k14=4.5
MaedaAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}

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
dt=0.001
t_tot = (30+30)*Nt_Maeda
t=np.arange(0,t_tot,dt)

# define extracellular cAMP input amounts 
cAMPe_in = np.array([0.005, 0.05, 0.5])   # np.logspace(-1, 2, num=6) # 

# preallocate trace arrays
cAMPi_traces = np.zeros([len(cAMPe_in),int(len(t)/2)])

# iterate through each eCAMP amount 
for i in range(len(cAMPe_in)):
    
    # define cAMPe input
    camp_input_trace=np.zeros(len(t))
    camp_input_trace[int(round(0.75*t_tot/dt)):] = cAMPe_in[i]
    
    # initialize population 
    pop = mlf.Population([1,1],initialVals,MaedaAgentParam,t)
    
    # run simulation 
    pop.run(dt,camp_input_trace)
    
    # store cAMPi trace value
    cAMPi_traces[i] = pop.cAMPi[int(len(t)/2):]/Nh_Maeda

# save adjusted time trace 
t = np.arange(0,t_tot/2,dt)/Nt_Maeda

# export time and cAMPi data
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/popAddcAMP_data/popAddcAMP_Maeda.pickle','wb') as f:
    pickle.dump([t,cAMPi_traces],f,pickle.HIGHEST_PROTOCOL)


#%% Kamino 2017, fig 5D group oscillations
import Kamino2017Func as kf

# set model parameters
tau=1.5
n=2
K=4
kt=2
delta=0.01
gamma = 3
rho = 1 # cell density
params={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
       'gamma':gamma,'rho':rho}

# set initial values
x0=0.01
y0=0.08
z0=0.01

# set simulation time parameters
dt=0.001
t_tot=30 *Nt_Kamino
t = np.arange(0,t_tot,dt)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              

# define cAMPe amounts  
cAMPe_in = np.array([0.001,0.002,0.1])

# preallocate trace arrays
x_traces = np.zeros([len(cAMPe_in),len(t)])
y_traces = np.zeros_like(x_traces)
z_traces = np.zeros_like(x_traces)

# iterate through each cAMPe amount and run population
for i in range(len(cAMPe_in)):
    
    # define cAMPe trace 
    camp_input_trace=np.zeros(len(t))
    camp_input_trace[int(round(0.5*t_tot/dt)):] = cAMPe_in[i]
    
    # initialize cell population 
    pop = kf.Population([x0,y0,z0],params,t)
    
    # run simulation
    pop.run(camp_input_trace,dt)

    # save traces
    x_traces[i] = pop.x
    y_traces[i] = pop.y/Nh_Kamino # normalize y trace 
    z_traces[i] = pop.z

# save adjusted t trace 
t = t/Nt_Kamino

# save to pickle file 
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/popAddcAMP_data/popAddcAMP_Kamino.pickle','wb') as f:
    pickle.dump([t,x_traces,y_traces,z_traces],f,pickle.HIGHEST_PROTOCOL)

