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

#%% Fig1 E F G & H
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
    
# Traces
A_trace_offset=1.5
A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
A_trace_plot=(A_trace_orig+A_trace_offset)/Na;

R_trace_orig = np.array(R_trace_orig) 
t_plot = np.array(t)/Nt
'''
# nullclines
A_null=np.arange(-2.5,2.5,0.025)
dAdt_null_no_stim=A_null-1/3*np.power(A_null,3)#  dA/dt=0 nullcline w/o stim
a=test_agent.AgentParam['e']
Kd= test_agent.AgentParam['Kd']
dAdt_null_stim=A_null-1/3*np.power(A_null,3)+a*np.log(1+cAMP/Kd) # dA/dt=0 nullcline with stim
c0= test_agent.AgentParam['c0']
g= test_agent.AgentParam['g'] # gamma
dRdt_null=1/g*(A_null+c0) # dR/dt=0 nullcline

fig,ax=plt.subplots(1,2)

ax[0].plot(t_plot,A_trace_plot)
ax[0].set_title('Activator after '+str(cAMP)+'nM cAMP stimulation')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Activator')

ax[1].plot(A_null,dRdt_null,'r')
ax[1].plot(A_null,dAdt_null_no_stim,'g')
ax[1].plot(A_null,dAdt_null_stim,'k')
ax[1].plot(A_trace_orig,R_trace_orig,'k')
ax[1].set_ylim([-1, 2.5])
ax[1].set_xlabel('Activator')
ax[1].set_ylabel('Repressor')
plt.show()
'''
# Plot adaptive spike
label_font_size = 25
trace_width = 6.0

fig5 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot,signal_trace, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
# ax1.set_xlabel('Time',fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot,A_trace_plot, color='g',linewidth=trace_width)
ax2.set_ylabel(r'Activator',fontsize=label_font_size)
ax3.yaxis.label.set_color('g')

ax3 = fig5.add_subplot(grid[1:, 0], sharex=ax2, frameon=False)
line2=ax3.plot(t_plot,R_trace_orig, color='b', linewidth=trace_width) # right axis
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")
ax3.set_ylabel('Repressor' ,fontsize=label_font_size)
ax3.yaxis.label.set_color('b')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax3.set_xlabel('Time',fontsize=label_font_size)

plt.show()


#%% Fig1 E & F- nullclines


