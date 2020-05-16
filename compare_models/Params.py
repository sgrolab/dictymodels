# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 11:00:54 2020

@author: Chuqiao Huyan

Time and cAMP maginitude normalization parameters. These parameters are kept
constant across the whole paper for each model.

"""
import math

# Normalization parameters for time scale and response height
NormParams = {'Nt_Gregor':6, 'Nt_Sgro':27, 'Nt_Goldbeter':6.94, 
             'Nt_Maeda':3.57, 'Nt_Kamino': 5.23,
             'Nh_Gregor': 19.6, 'Nh_Sgro': 3.5, 'Nh_Goldbeter':210.53,
             'Nh_Maeda':3.15, 'Nh_Kamino':0.26,
             'Nh_Sgro_offset':-1.5, 'Nh_Kamino_offset': 0.0585}


# Martiel 1987, Table 2 parameters 
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

# Goldbeter3PopParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
#             'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
#             'ki':ki,'kt':kt, 'kc':kc,'h':h}
sigma = 10 # noise strength
N = 100 # number of cells
Goldbeter3PopParam = {'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h, 'sigma':sigma, 'N':N}

############
## Laub loomis 1998 parameters
#k1=1.4; k2=0.9; k3=2.5; k4=1.5; k5=0.6
#k6=0.8; k7=2.0; k8=1.3; k9=0.7; k10=1.0
#k11=0.3; k12=3.1; k13=1.8; k14=1.5

# Maeda & Loomis 2004 parameters
# parameters from Maeda & Loomis 2004 paper
k1=2.0; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=1.0; k8=1.3; k9=0.3; k10=0.8
k11=0.7; k12=4.9; k13=23; k14=4.5
MaedaAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}
    
sigma = 0.1 # noise strength
N = 100  # number of cells
MaedaPopParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14, 'N':N, 'sigma':sigma}

# Gregor 2010 parameters
Amax=20;  Abas=0.4 # uM
w=2*math.pi/6 # min-1
Vc=1.1e-9 # ml
St=1.33 # cm2
Sc=1.3e-6 # cm2
K=0.0004 # uM, 400 pM
c_sec= 3.6 # min-1
c_excite=1.01 # min-1
eta=0.02 # noise stength
GregorAgentParam={'Amax':Amax,'Abas':Abas,'w':w,'Vc':Vc,'St':St,'Sc':Sc,'K':K,\
            'c_sec':c_sec,'c_excite':c_excite,'eta':eta}

Nc=100 # Num of cells
rho = 1/12
k = 5 #ml/min
Vt = 1 #chamber size ml
GregorPopParam={'Amax':Amax,'Abas':Abas,'w':w,'Vc':Vc,'St':St,'Sc':Sc,'K':K,\
            'c_sec':c_sec,'c_excite':c_excite,'Vt':Vt,'Nc':Nc,'eta':eta,\
                'rho':rho,'k':k}
    
# Sgro 2015 parameters
SgroAgentParam={'e':0.1,'tauA':0.09,'tauR':0.9,'g':0.5,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
    
e=0.1; tauA=0.09; tauR=tauA/e; g=0.5; sigma = 0.15 # noise strength
N = 100 # number of cells in the population
rho = 10**(-3.5); j = 0.5
SgroPopParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':sigma,'N':N,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'offset_A':1.5,'flux_thrs':0, 'rho': rho,'j': j}

    
# Kamino 2017 parameters
tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
KaminoAgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
    
sigma = 0.01# noise strength
N=100 # number of cells in a population
KaminoPopParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,'sigma':sigma, 'N':N,'rho':rho,'gamma':gamma}