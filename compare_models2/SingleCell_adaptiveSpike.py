# Single Cell, Adaptive Spike 

import os, sys, pickle 
import numpy as np 
import matplotlib.pyplot as plt 

# add path for function files
sys.path.insert(0,'//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

# import normalization parameters
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/normParams.pickle','rb') as f:
    Nt_Gregor,Nt_Sgro,Nt_Goldbeter,Nt_Maeda,Nt_Kamino,Nh_Gregor,Nh_Sgro,Nh_Goldbeter,Nh_Maeda,Nh_Kamino,Nh_Kamino_offset,Nh_Sgro_offset = pickle.load(f)


#%% Gregor 

import Gregor2010Func as gf

# set parameters
Amax=20
Abas=0.4            # uM
w = np.pi/3       # min-1
Vc=1.1e-9           # ml
St=1.33             # cm2
Sc=1.3e-6           # cm2
K=0.0004            # uM, 400 pM
c_sec= 3.6          # min-1
c_excite=1.01       # min-1
GregorAgentParam={'Amax':Amax,'Abas':Abas,'w':w,'Vc':Vc,'St':St,'Sc':Sc,'K':K,\
            'c_sec':c_sec,'c_excite':c_excite}

# set simulation time parameters
dt=0.001
t_tot=6*Nt_Gregor
t=np.arange(0,t_tot,dt)

# compute initial values 
cAMPi0 = 0.4
thetai0 = np.arcsin((cAMPi0*2-Amax-Abas)/(-Amax+Abas))
cAMPe0 = 0

# define cAMPe_in trace
constant_signal=0
signal_trace=np.zeros(len(t))
signal_trace[int(len(t)/6):] = constant_signal

# set noise amount 
eta = 0 # noise stength
rng = np.random.default_rng(seed=1)
r = np.sqrt(dt) * rng.uniform(-1,1,size = (len(t))) * eta

# initialize agent
cell = gf.Cell([1,1],[cAMPi0,thetai0,cAMPe0],GregorAgentParam,t)

# run simulation 
cell.run(dt,r,signal_trace) 

plt.plot(cell.t/Nt_Gregor,cell.cAMPi/Nh_Gregor)
plt.ylim([-0.5,1.2])

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_adaptiveSpike_data/adaptiveSpike_Gregor.pickle','wb') as f:
    pickle.dump([cell.t/Nt_Gregor,cell.cAMPi/Nh_Gregor],f,pickle.HIGHEST_PROTOCOL)


#%% Sgro 

import Sgro2015Func as sf

# set parameters 
SgroAgentParam={'e':0.1,'tauA':0.09,'tauR':0.9,'g':0.5,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}

# set time parameters
dt=0.005
t_tot=6*Nt_Sgro
t=np.arange(0,t_tot,dt)

# set initial values 
A0=-1.5
R0=-0.5

# define cAMPe_in trace 
constant_signal= 1
signal_trace=np.zeros(len(t))
signal_trace[int(len(t)/6):] = constant_signal

# set noise amount 
rng = np.random.default_rng(seed=1)
r = np.sqrt(dt) * rng.normal(0,1,size = (len(t)))

# initialize cell 
cell = sf.Cell([1,1],[A0,R0],SgroAgentParam,t)

# run simulation 
cell.run(dt,signal_trace,r)

plt.plot(cell.t/Nt_Sgro,(cell.A-Nh_Sgro_offset)/Nh_Sgro)
plt.ylim([-0.5,1.2])
plt.xlim([0,6])

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_adaptiveSpike_data/adaptiveSpike_Sgro.pickle','wb') as f:
    pickle.dump([cell.t/Nt_Sgro,(cell.A-Nh_Sgro_offset)/Nh_Sgro],f,pickle.HIGHEST_PROTOCOL)


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
t_tot=6*Nt_Goldbeter
t=np.arange(0,t_tot,dt)
    
# set initial values 
p0=0.8
a0=3
b0=0.9
g0=0

# define cAMPe_in trace 
constant_signal= 1
signal_trace=np.zeros(len(t))
signal_trace[int(len(t)/6):] = constant_signal

# initialize cell 
cell = gbf.Cell(Goldbeter3AgentParam,[p0,b0,g0],t)

# run cell 
cell.run(dt,a0,signal_trace)

plt.plot(cell.t/Nt_Goldbeter,cell.b/Nh_Goldbeter)
plt.ylim([-0.5,1.2])
plt.xlim([0,6])

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_adaptiveSpike_data/adaptiveSpike_Goldbeter.pickle','wb') as f:
    pickle.dump([cell.t/Nt_Goldbeter,cell.b/Nh_Goldbeter],f,pickle.HIGHEST_PROTOCOL)


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
dt=0.0001
t_tot=6*Nt_Maeda
t=np.arange(0,t_tot,dt)

# set initial values 
ACA0=0.1
PKA0=0.1
ERK20=0.1
RegA0=0.1
cAMPi0=0.01
cAMPe0=0.1
CAR10=0.1

# define cAMPe_in trace 
constant_signal = 1
signal_trace=np.zeros(len(t))
signal_trace[int(len(t)/6):] = constant_signal

# initialize cell 
cell = mlf.Cell([0,0],[ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10],MaedaAgentParam,t)

# run cell
cell.run(dt,signal_trace)

plt.plot(cell.t/Nt_Maeda,cell.cAMPi/Nh_Maeda)
plt.xlim([0,6])
plt.ylim([-0.5,1.2])

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_adaptiveSpike_data/adaptiveSpike_MaedaLoomis.pickle','wb') as f:
    pickle.dump([cell.t/Nt_Maeda,cell.cAMPi/Nh_Maeda],f,pickle.HIGHEST_PROTOCOL)


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
dt=0.0001
t_tot=6*Nt_Kamino
t = np.arange(0,t_tot,dt)

# set initial values 
x0=0.01
y0=0.06
z0=0.005

# define cAMPe_in trace 
constant_signal = 1 
signal_trace=np.zeros(len(t))
signal_trace[int(len(t)/6):] = constant_signal

# initialize cell 
cell = kf.Cell([x0,y0,z0],KaminoAgentParam,t)

# run cell
cell.run(dt,signal_trace)

plt.plot(cell.t/Nt_Kamino,(cell.y-Nh_Kamino_offset)/Nh_Kamino)
plt.xlim([0,6])
plt.ylim([-0.5,1.2])

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_adaptiveSpike_data/adaptiveSpike_Kamino.pickle','wb') as f:
    pickle.dump([cell.t/Nt_Kamino,(cell.y-Nh_Kamino_offset)/Nh_Kamino],f,pickle.HIGHEST_PROTOCOL)

