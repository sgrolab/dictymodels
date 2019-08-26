# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 11:53:44 2019

@author: ellin

Maeda and Loomis model from 2004Sience paper

"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt
# remove all variables in the work space
# sys.modules[__name__].__dict__.clear()

class LaubLoomis1998_agent:
    def __init__(self,pos,state,AgentParam):
        self.pos=pos
        self.state=state
        self.AgentParam=AgentParam
        self.ACA_now=state[0]
        self.PKA_now=state[1]
        self.ERK2_now=state[2]
        self.RegA_now=state[3]
        self.cAMPi_now=state[4]
        self.cAMPe_now=state[5]
        self.CAR1_now=state[6] # initial state input as a list variable [ACA0,PKA0,....]
        
    
    def update(self,dt,overriding_sig='none'):
        k1=AgentParam['k1']   
        k2=AgentParam['k2']  
        k3=AgentParam['k3']  
        k4=AgentParam['k4']  
        k5=AgentParam['k5']  
        k6=AgentParam['k6']  
        k7=AgentParam['k7']  
        k8=AgentParam['k8']  
        k9=AgentParam['k9']  
        k10=AgentParam['k10']  
        k11=AgentParam['k11']  
        k12=AgentParam['k12']  
        k13=AgentParam['k13']  
        k14=AgentParam['k14']  
     
        
        ACA_next=self.ACA_now+(k1*self.CAR1_now-k2*self.ACA_now*self.PKA_now)*dt
        PKA_next=self.PKA_now+(k3*self.cAMPi_now-k4*self.PKA_now)*dt
        ERK2_next=self.ERK2_now+(k5*self.CAR1_now-k6*self.PKA_now*self.ERK2_now)*dt
        RegA_next=self.RegA_now+(k7-k8*self.ERK2_now*self.RegA_now)*dt
        cAMPi_next=self.cAMPi_now+(k9*self.ACA_now-k10*self.RegA_now*self.cAMPi_now)*dt
        
        if isinstance(overriding_sig, str):
            cAMPe_next=self.cAMPe_now+(k11*self.ACA_now-k12*self.cAMPe_now)*dt
        else:
            cAMPe_next=overriding_sig
            
        CAR1_next=self.CAR1_now+(k13*self.cAMPe_now-k14*self.CAR1_now)*dt
        
        
        self.ACA_now=ACA_next
        self.PKA_now=PKA_next
        self.ERK2_now=ERK2_next
        self.RegA_now=RegA_next
        self.cAMPi_now=cAMPi_next
        self.cAMPe_now=cAMPe_next
        self.CAR1_now=CAR1_next
        # update the current states
        
        return ACA_next,PKA_next,ERK2_next,RegA_next,cAMPi_next,cAMPe_next,CAR1_next
    """
    def flux(self,signals):
        if self.A>self.AgentParam['flux_thrs']:
            agent_flux=True
        else:
            agent_flux=False
        return agent_flux
    """

    def print_state(self):
        print('past A:'+str(self.A_now))
        print('past R:'+str(self.R_now)) 
        
#%%
k1=1.4
k2=0.9
k3=2.5
k4=1.5
k5=0.6
k6=0.8
k7=2.0
k8=1.3
k9=0.7
k10=1.0
k11=0.3
k12=3.1
k13=1.8
k14=1.5
AgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}
ACA0=0.1
PKA0=0.1
ERK20=0.1
RegA0=0.1
cAMPi0=0.01
cAMPe0=0.1
CAR10=0.1


state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]
LaubLoomis_agent=LaubLoomis1998_agent([1,1],state0,AgentParam)

# Simulate time traces

"""
# defining signals
cAMP=1
tTrig=int(round(1/dt*Nt)) # Time step at which triggering happens
signal_trace=[0]*len(t)
signal_trace[tTrig:]=[cAMP]*len(signal_trace[tTrig:])

# test update function with one step caculations
# A,R,r=test_agent.update(1,dt)
"""
# initializations
ACA_trace=[ACA0]
PKA_trace=[PKA0]
ERK2_trace=[ERK20]
RegA_trace=[RegA0]
cAMPi_trace=[cAMPi0]
# cAMPe_trace=[cAMPe0]
CAR1_trace=[CAR10]

dt=0.001
t_tot=25
t=list(np.arange(0,t_tot,dt))

stim_time_step=3000 # at this time step input is applied
constant_signal=5; # randomly set 2019/7/17
signal_trace=np.zeros(len(t))
signal_trace[stim_time_step:] = constant_signal

for i in range(len(t)-1):
    ACA_now=ACA_trace[i]
    PKA_now=PKA_trace[i]
    ERK2_now=ERK2_trace[i]
    RegA_now=RegA_trace[i]
    cAMPi_now=cAMPi_trace[i]
    cAMPe_now=cAMPi_trace[i]
    CAR1_now=CAR1_trace[i]
    
    ACA_next,PKA_next,ERK2_next,RegA_next,\
    cAMPi_next,cAMPe_next,CAR1_next=LaubLoomis_agent.update(dt,signal_trace[i])
    
    ACA_trace.append(ACA_next)
    PKA_trace.append(PKA_next)
    ERK2_trace.append(ERK2_next)
    RegA_trace.append(RegA_next)
    cAMPi_trace.append(cAMPi_next)
    # cAMPe_trace.append(cAMPe_next)
    CAR1_trace.append(CAR1_next)
    

ERK2_trace = np.array(ERK2_trace) # vectorize p_trace
cAMPi_trace = np.array(cAMPi_trace)
t_plot = np.array(t)


# Figure single cell adaptive spike 
label_font_size = 25
trace_width = 6.0

fig5 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot,signal_trace, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
# ax1.set_xlabel('Time',fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot,cAMPi_trace, color='g',linewidth=trace_width)
ax2.set_ylabel(r'$cAMP_{cyto}$',fontsize=label_font_size)
ax2.yaxis.label.set_color('g')

ax3 = fig5.add_subplot(grid[1:, 0], sharex=ax2, frameon=False)
line2=ax3.plot(t_plot,ERK2_trace, color='b', linewidth=trace_width) # right axis
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")
ax3.set_ylabel('ERK2' ,fontsize=label_font_size)
ax3.yaxis.label.set_color('b')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax3.set_xlabel('Time',fontsize=label_font_size)

plt.show()


