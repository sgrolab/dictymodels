# Single Cell, Fold Change Detection 

import sys, pickle 
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.signal import find_peaks

# add path for function files
sys.path.insert(0,'//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

# import normalization parameters
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/normParams.pickle','rb') as f:
    Nt_Gregor,Nt_Sgro,Nt_Goldbeter,Nt_Maeda,Nt_Kamino,Nh_Gregor,Nh_Sgro,Nh_Goldbeter,Nh_Maeda,Nh_Kamino,Nh_Kamino_offset,Nh_Sgro_offset = pickle.load(f)

#%% Kamino 2017 (IFFL)

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
dt = 0.001
t_tot = 30 * Nt_Kamino
t = np.arange(0,t_tot,dt)

# set initial values 
x0=0.01
y0=0.06
z0=0.005

# define cAMPe amounts 
cAMPe_in_priming = np.array((0.1,1,3,10))
cAMPe_in_foldChange = np.logspace(0,2,8)

# preallocate second peak prominence array 
peak2Prom = np.zeros([len(cAMPe_in_priming),len(cAMPe_in_foldChange)])

# iterate through each priming concentration and fold change 
for i in range(len(cAMPe_in_priming)):
    for j in range(len(cAMPe_in_foldChange)):
        
        # define cAMPe_in trace 
        cAMPe_in = np.zeros(len(t))
        cAMPe_in[int(len(t)/3):] = cAMPe_in_priming[i]
        cAMPe_in[int(2*len(t)/3):] = cAMPe_in_priming[i] * cAMPe_in_foldChange[j]
        
        # initialize agent 
        cell = kf.Cell([x0,y0,z0],KaminoAgentParam,t)
        
        # run simulation 
        cell.run(dt,cAMPe_in)

        # get cAMPi traces during cAMPe exposures
        cAMPi_1 = (cell.y[int(len(t)/3):int(2*len(t)/3)]-Nh_Kamino_offset)/Nh_Kamino
        cAMPi_2 = (cell.y[int(2*len(t)/3):]-Nh_Kamino_offset)/Nh_Kamino
        
        # detect peaks 
        Pk1Pos, Pk1Props = find_peaks(cAMPi_1, prominence=(0,50/Nh_Kamino))
        Pk2Pos, Pk2Props = find_peaks(cAMPi_2, prominence=(0,50/Nh_Kamino))
        
        # if second spike is present, save values 
        if Pk2Pos.size > 0:
            peak2Prom[i,j]=Pk2Props["prominences"][0]

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_foldChange_data/FCD_Kamino.pickle','wb') as f:
    pickle.dump(peak2Prom,f,pickle.HIGHEST_PROTOCOL)

#%% Sgro 2015 (IPNFB) 

import Sgro2015Func as sf

# set parameters 
SgroAgentParam={'e':0.1,'tauA':0.09,'tauR':0.9,'g':0.5,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}

# set time parameters
dt=0.005
t_tot=30*Nt_Sgro
t=np.arange(0,t_tot,dt)

# set initial values 
A0=-1.5
R0=-0.5

# initialize random number generator 
rng = np.random.default_rng(seed=1)

# define cAMPe amounts 
cAMPe_in_priming = np.array((0.1,1,3,10))
cAMPe_in_foldChange = np.logspace(0,2,8)

# set number of reptitions 
nCells = 50

# preallocate array 
peak2Prom = np.zeros([len(cAMPe_in_priming),len(cAMPe_in_foldChange),nCells])

# iterate through each priming concentration and fold change 
for i in range(len(cAMPe_in_priming)):
    for j in range(len(cAMPe_in_foldChange)):
        for k in range(nCells):
            
            # define cAMPe_in trace 
            cAMPe_in = np.zeros(len(t))
            cAMPe_in[int(len(t)/3):] = cAMPe_in_priming[i]
            cAMPe_in[int(2*len(t)/3):] = cAMPe_in_priming[i] * cAMPe_in_foldChange[j]
            
            # initialize cell 
            cell = sf.Cell([1,1],[A0,R0],SgroAgentParam,t)
            
            # generate noise trace 
            r = np.sqrt(dt) * rng.normal(0,1,size = (len(t)))
            
            # run simulation 
            cell.run(dt,cAMPe_in,r)
            
            # get cAMPi traces for one spike after cAMPe change 
            cAMPi_1 = (cell.A[int(len(t)/3):int(len(t)/3+1.5*Nt_Sgro/dt)]-Nh_Sgro_offset)/Nh_Sgro
            cAMPi_2 = (cell.A[int(2*len(t)/3):int(2*len(t)/3+1.5*Nt_Sgro/dt)]-Nh_Sgro_offset)/Nh_Sgro
            
            # detect peaks 
            Pk1Pos, Pk1Props = find_peaks(cAMPi_1, prominence=(0.4,2))
            Pk2Pos, Pk2Props = find_peaks(cAMPi_2, prominence=(0.4,2))
            
            # if second spike is present, save values 
            if Pk2Pos.size > 0:
                peak2Prom[i,j,k]=Pk2Props["prominences"][0]

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_foldChange_data/FCD_Sgro.pickle','wb') as f:
    pickle.dump(peak2Prom,f,pickle.HIGHEST_PROTOCOL)

#%% Goldbeter (Receptor Desensitization) 

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
dt=0.0005
t_tot=30*Nt_Goldbeter
t=np.arange(0,t_tot,dt)

# set initial values 
p0=0.8
a0=3
b0=0.9
g0=0

# define cAMPe amounts 
cAMPe_in_priming = np.array((0.1,1,3,10))
cAMPe_in_foldChange = np.logspace(0,2,12)

# preallocate second peak prominence array 
peak2Prom = np.zeros([len(cAMPe_in_priming),len(cAMPe_in_foldChange)])

# iterate through each priming concentration and fold change 
for i in range(len(cAMPe_in_priming)):
    for j in range(len(cAMPe_in_foldChange)):
        
        # define cAMPe_in trace 
        cAMPe_in = np.zeros(len(t))
        cAMPe_in[int(len(t)/3):] = cAMPe_in_priming[i]
        cAMPe_in[int(2*len(t)/3):] = cAMPe_in_priming[i] * cAMPe_in_foldChange[j]
        
        # initialize cell 
        cell = gbf.Cell(Goldbeter3AgentParam,[p0,b0,g0],t)
        
        # run cell 
        cell.run(dt,a0,cAMPe_in)
        
        # get cAMPi traces after cAMPe changes
        cAMPi_1 = cell.b[int(len(t)/3):int(2*len(t)/3)]/Nh_Goldbeter
        cAMPi_2 = cell.b[int(2*len(t)/3):]/Nh_Goldbeter
        
        # detect peaks 
        Pk1Pos, Pk1Props = find_peaks(cAMPi_1, prominence=(0,3000/Nh_Goldbeter))
        Pk2Pos, Pk2Props = find_peaks(cAMPi_2, prominence=(0,3000/Nh_Goldbeter))
        
        # if second spike is present, save values 
        if Pk2Pos.size > 0:
            peak2Prom[i,j]=Pk2Props["prominences"][0]

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_foldChange_data/FCD_Goldbeter.pickle','wb') as f:
    pickle.dump(peak2Prom,f,pickle.HIGHEST_PROTOCOL)

#%% sample for Fig 6A

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
dt = 0.001
t_tot = 30 * Nt_Kamino
t = np.arange(0,t_tot,dt)

# set initial values 
x0=0.01
y0=0.06
z0=0.005

# define cAMPe amounts 
cAMPe_in_priming = 1
cAMPe_in_foldChange = 5


# define cAMPe_in trace 
cAMPe_in = np.zeros(len(t))
cAMPe_in[int(len(t)/3):] = cAMPe_in_priming
cAMPe_in[int(2*len(t)/3):] = cAMPe_in_priming * cAMPe_in_foldChange

# initialize agent 
cell = kf.Cell([x0,y0,z0],KaminoAgentParam,t)

# run simulation 
cell.run(dt,cAMPe_in)

# save trace 
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_foldChange_data/FCD_sample.pickle','wb') as f:
    pickle.dump([cell.t,cell.y],f,pickle.HIGHEST_PROTOCOL)