# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:31:57 2019

@author: Chuqiao
"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt
from Sgro2015_agent import *
# import sys
# remove all variables in the work space
# sys.modules[__name__].__dict__.clear() 
e=0.1
tauA=0.09
tauR=tauA/e
g=0.5
AgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
Nt=27; # normalization factor of t
Na=3.5;  # normalization factor of A
A0=-1.5
R0=-0.5
test_agent=Sgro2015_agent([1,1],[A0,R0],AgentParam)

#%% Fig1 G & H
dt=0.005 
t_tot=175
t=list(np.arange(0,t_tot,dt))

# defining signals
cAMP=1
tTrig=int(round(1/dt*Nt)) # Time step at which triggering happens
signal_trace=[0]*len(t)
signal_trace[tTrig:]=[cAMP]*len(signal_trace[tTrig:])

# test update function with one step caculations
# A,R,r=test_agent.update(1,dt)

# initializations
A_trace_orig=[A0]
R_trace_orig=[R0]
r_trace=[]

for i in range(len(t)-1):
    A_now=A_trace_orig[i]
    R_now=R_trace_orig[i]
    signal_now=signal_trace[i]
    
    A_next,R_next,r_now=test_agent.update(signal_now,dt)
    A_trace_orig.append(A_next)
    R_trace_orig.append(R_next)
    r_trace.append(r_now)
    

A_trace_offset=1.5
A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
A_trace_plot=(A_trace_orig+A_trace_offset)/Na;
t_plot = np.array(t)/Nt

plt.plot(t_plot,A_trace_plot)
plt.xlabel('Time')
plt.ylabel('Activator')
plt.show()

#%% Fig1 E & F
A_null=np.arange(-2.5,2.5,0.025)

