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

class Maeda2004_agent:
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
        
    
    def update(self,signals,dt):
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
        cAMPe_next=self.cAMPe_now+(k11*self.ACA_now-k12*self.cAMPe_now)*dt
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
        

k1=2.0
k2=0.9
k3=2.5
k4=1.5
k5=0.6
k6=0.8
k7=1.0
k8=1.3
k9=0.3
k10=0.8
k11=0.7
k12=4.9
k13=23.0
k14=4.5
AgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}
ACA0=0
PKA0=0
ERK20=0.4
RegA0=0
cAMPi0=0.5
cAMPe0=0.1
CAR10=0

[0,0,0.4,0,0.5,0.1,0] #
state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]
test_agent=Maeda2004_agent([1,1],state0,AgentParam)

#%% Simulate time traces
dt=0.01  
t_tot=40
t=list(np.arange(0,t_tot,dt))
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
ACA=[ACA0]
PKA=[PKA0]
ERK2=[ERK20]
RegA=[RegA0]
cAMPi=[cAMPi0]
cAMPe=[cAMPe0]
CAR1=[CAR10]


for i in range(len(t)-1):
    ACA_now=ACA[i]
    PKA_now=PKA[i]
    ERK2_now=ERK2[i]
    RegA_now=RegA[i]
    cAMPi_now=cAMPi[i]
    cAMPe_now=cAMPi[i]
    CAR1_now=CAR1[i]
    
    signal_now=0
    
    ACA_next,PKA_next,ERK2_next,RegA_next,\
    cAMPi_next,cAMPe_next,CAR1_next=test_agent.update(signal_now,dt)
    
    ACA.append(ACA_next)
    PKA.append(PKA_next)
    ERK2.append(ERK2_next)
    RegA.append(RegA_next)
    cAMPi.append(cAMPi_next)
    cAMPe.append(cAMPe_next)
    CAR1.append(CAR1_next)
    

plt.plot(t,cAMPe)
plt.hold(True)
plt.plot(t,ERK2)
plt.hold(True)
plt.plot(t,cAMPi)
plt.xlabel('Time (min)')
plt.ylabel('Stuff')
plt.show()


