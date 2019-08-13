# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 16:59:30 2019

Sgro 2015 population model, single cells as agents

@author: ellin
"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal
from scipy.signal import find_peaks

# remove all variables in the work space
# sys.modules[__name__].__dict__.clear()

def update(Param,state,signals,dt):
    e=Param['e']   
    tauA=Param['tauA'] 
    tauR=Param['tauR']
    g=Param['g'] 
    c0=Param['c0'] 
    sigma=Param['sigma']
    a=Param['a'] 
    alpha0=Param['alpha0'] 
    alpha_pde=Param['alpha_pde'] 
    Kd=Param['Kd'] 
    S=Param['S'] 
    
    A_now=state[0]
    R_now=state[1]
    
    fA=A_now-(A_now)**3/3-R_now+a*np.log(1+signals/Kd)
    fR=e*(A_now-g*R_now+c0)
    r= math.sqrt(dt)*random.gauss(0, 1)
    A_next=A_now+fA*dt+sigma*r
    R_next=R_now+fR*dt/tauR
    return A_next,R_next,r
#%%
NumofCells=100;

A_pop_0=(-3)*np.random.rand(NumofCells,1)  
R_pop_0=(-1)*np.random.rand(NumofCells,1)  

# single cell parameters
e=0.1
tauA=0.09
tauR=tauA/e
g=0.5
Param={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
Nt=27; # normalization factor of t
Na=3.5;  # normalization factor of A

alphaf=2.5
rho=10**(-3.5)
J=0.6

dt=0.005 
t_tot=175
t=list(np.arange(0,t_tot,dt))

cAMP_ex=np.zeros((1,len(t)))

A_pop = np.zeros((NumofCells,len(t))) 
A_pop[:,:1]=A_pop_0
R_pop = np.zeros((NumofCells,len(t))) 
R_pop[:,:1]=R_pop_0
r_pop=np.zeros((NumofCells,len(t))) 
'''
pop = list() # store population as a list of agents
for i in range(NumofCells):
    A0=A_pop_0[i]
    R0=R_pop_0[i]
    pop.append(Sgro2015_agent([1,1],[A0,R0],AgentParam))
'''

for i in range(len(t)-1):
    cAMP_ex_now=cAMP_ex[0][i]
    for j in range(NumofCells):
        
        A_now=A_pop[j][i]
        R_now=R_pop[j][i]
        
        A_next,R_next,r_now=update(Param,[A_now, R_now], cAMP_ex_now,dt)
        A_pop[j][i+1]=A_next
        R_pop[j][i+1]=R_next
        r_pop[j][i]=r_now
    
    plt.show() # set up figure out side of loop, plots in the loop like a video
    if (np.mod(i,10000) == 0): # check activator level every 500 time points
        A_mesh=np.reshape(A_pop[:,i],(10,10))
        plt.pcolor(A_mesh)
        plt.clim(-2.5,2.5)
        plt.title('Activator value in the population at time: %i' %(i*dt))
        plt.pause(0.02)
    
############### here 4/25 
    
#%%


