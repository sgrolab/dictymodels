# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 10:10:49 2019

@author: Chuqiao
"""

import numpy as np
import random
import math
import matplotlib.pyplot as plt
# remove all variables in the work space
# sys.modules[__name__].__dict__.clear()

class MaedaLoomis2004_agent:
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
        
    
    def update(self,dt,signal):
        k1=self.AgentParam['k1']   
        k2=self.AgentParam['k2']  
        k3=self.AgentParam['k3']  
        k4=self.AgentParam['k4']  
        k5=self.AgentParam['k5']  
        k6=self.AgentParam['k6']  
        k7=self.AgentParam['k7']  
        k8=self.AgentParam['k8']  
        k9=self.AgentParam['k9']  
        k10=self.AgentParam['k10']  
#        k11=self.AgentParam['k11']  
#        k12=self.AgentParam['k12']  
        k13=self.AgentParam['k13']  
        k14=self.AgentParam['k14']  
     
        
        ACA_next=self.ACA_now+(k1*self.CAR1_now-k2*self.ACA_now*self.PKA_now)*dt
        PKA_next=self.PKA_now+(k3*self.cAMPi_now-k4*self.PKA_now)*dt
        ERK2_next=self.ERK2_now+(k5*self.CAR1_now-k6*self.PKA_now*self.ERK2_now)*dt
        RegA_next=self.RegA_now+(k7-k8*self.ERK2_now*self.RegA_now)*dt
        cAMPi_next=self.cAMPi_now+(k9*self.ACA_now-k10*self.RegA_now*self.cAMPi_now)*dt
        cAMPe_next=signal
            
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
        pass

class MaedaLoomis2004_pop:
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
        
    
    def update(self,dt,campExt_influx,rho = 1,gamma = 0):
        k1=self.AgentParam['k1']   
        k2=self.AgentParam['k2']  
        k3=self.AgentParam['k3']  
        k4=self.AgentParam['k4']  
        k5=self.AgentParam['k5']  
        k6=self.AgentParam['k6']  
        k7=self.AgentParam['k7']  
        k8=self.AgentParam['k8']  
        k9=self.AgentParam['k9']  
        k10=self.AgentParam['k10']  
        k11=self.AgentParam['k11']  
        k12=self.AgentParam['k12']  
        k13=self.AgentParam['k13']  
        k14=self.AgentParam['k14']  
     
        
        ACA_next=self.ACA_now+(k1*self.CAR1_now-k2*self.ACA_now*self.PKA_now)*dt
        PKA_next=self.PKA_now+(k3*self.cAMPi_now-k4*self.PKA_now)*dt
        ERK2_next=self.ERK2_now+(k5*self.CAR1_now-k6*self.PKA_now*self.ERK2_now)*dt
        RegA_next=self.RegA_now+(k7-k8*self.ERK2_now*self.RegA_now)*dt
        cAMPi_next=self.cAMPi_now+(k9*self.ACA_now-k10*self.RegA_now*self.cAMPi_now)*dt
        
        cAMPe_next=self.cAMPe_now+(k11*rho*self.ACA_now-(k12+gamma)*(self.cAMPe_now - campExt_influx))*dt
        if cAMPe_next < 0 :
            cAMPe_next = 0
        CAR1_next=self.CAR1_now+(k13*self.cAMPe_now - k14*self.CAR1_now)*dt
        
        
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
        pass
    
#class LaubLoomis1998_pop:
#    def __init__(self,pos,state,AgentParam):
#        self.pos=pos
#        self.state=state
#        self.AgentParam=AgentParam
#        self.ACA_now=state[0]
#        self.PKA_now=state[1]
#        self.ERK2_now=state[2]
#        self.RegA_now=state[3]
#        self.cAMPi_now=state[4]
#        self.cAMPe_now=state[5]
#        self.CAR1_now=state[6] # initial state input as a list variable [ACA0,PKA0,....]
#        
#    
#    def update(self,dt,campExt_influx):
#        k1=self.AgentParam['k1']   
#        k2=self.AgentParam['k2']  
#        k3=self.AgentParam['k3']  
#        k4=self.AgentParam['k4']  
#        k5=self.AgentParam['k5']  
#        k6=self.AgentParam['k6']  
#        k7=self.AgentParam['k7']  
#        k8=self.AgentParam['k8']  
#        k9=self.AgentParam['k9']  
#        k10=self.AgentParam['k10']  
#        k11=self.AgentParam['k11']  
#        k12=self.AgentParam['k12']  
#        k13=self.AgentParam['k13']  
#        k14=self.AgentParam['k14']  
#     
#        
#        ACA_next=self.ACA_now+(k1*self.ERK2_now-k2*self.ACA_now)*dt
#        PKA_next=self.PKA_now+(k3*self.cAMPi_now-k4*self.PKA_now)*dt
#        ERK2_next=self.ERK2_now+(k5*self.CAR1_now-k6*self.PKA_now*self.ERK2_now)*dt
#        RegA_next=self.RegA_now+(k7-k8*self.ERK2_now*self.RegA_now)*dt
#        cAMPi_next=self.cAMPi_now+(k9*self.ACA_now-k10*self.RegA_now*self.cAMPi_now)*dt
#        
#        cAMPe_next=self.cAMPe_now+(k11*self.ACA_now-k12*(self.cAMPe_now - campExt_influx))*dt
#            
#        CAR1_next=self.CAR1_now+(k13*self.cAMPe_now-k14*self.CAR1_now*self.PKA_now)*dt
#        
#        
#        self.ACA_now=ACA_next
#        self.PKA_now=PKA_next
#        self.ERK2_now=ERK2_next
#        self.RegA_now=RegA_next
#        self.cAMPi_now=cAMPi_next
#        self.cAMPe_now=cAMPe_next
#        self.CAR1_now=CAR1_next
#        # update the current states
#        
#        return ACA_next,PKA_next,ERK2_next,RegA_next,cAMPi_next,cAMPe_next,CAR1_next
#    """
#    def flux(self,signals):
#        if self.A>self.AgentParam['flux_thrs']:
#            agent_flux=True
#        else:
#            agent_flux=False
#        return agent_flux
#    """
#
#    def print_state(self):
#        pass