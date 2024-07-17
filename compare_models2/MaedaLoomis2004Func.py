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

class Cell:
    def __init__(self,pos,initialVals,AgentParam,t):
        self.pos=pos
        self.AgentParam=AgentParam
        self.t = t 
        
        # preallocate variable arrays 
        self.ACA = np.zeros(len(t))
        self.PKA = np.zeros(len(t))
        self.ERK2 = np.zeros(len(t))
        self.RegA = np.zeros(len(t))
        self.cAMPi = np.zeros(len(t))
        self.cAMPe = np.zeros(len(t))
        self.CAR1 = np.zeros(len(t))
        
        # set initial values
        self.ACA[0] = initialVals[0]
        self.PKA[0] = initialVals[1]
        self.ERK2[0] = initialVals[2]
        self.RegA[0] = initialVals[3]
        self.cAMPi[0] = initialVals[4]
        self.cAMPe[0] = initialVals[5]
        self.CAR1[0] = initialVals[6] 
    
    def run(self,dt,cAMPe_in):
        
        # get parameter values from dictionary
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
        k13=self.AgentParam['k13']  
        k14=self.AgentParam['k14']
        
        # run simulation 
        for i in range(1,len(self.t)):
        
            # get variable values from previous timestep 
            ACA = self.ACA[i-1]
            PKA = self.PKA[i-1]
            ERK2 = self.ERK2[i-1]
            RegA = self.RegA[i-1]
            cAMPi = self.cAMPi[i-1]
            cAMPe = self.cAMPe[i-1]
            CAR1 = self.CAR1[i-1]
            
            # compute next values
            self.ACA[i] = ACA + (k1*CAR1 - k2*ACA*PKA)*dt
            self.PKA[i] = PKA + (k3*cAMPi - k4*PKA)*dt
            self.ERK2[i] = ERK2 + (k5*CAR1-k6*PKA*ERK2)*dt
            self.RegA[i] = RegA + (k7-k8*ERK2*RegA)*dt
            self.cAMPi[i] = cAMPi + (k9*ACA-k10*RegA*cAMPi)*dt
            self.cAMPe[i] = cAMPe_in[i]
            self.CAR1[i] = CAR1 + (k13*cAMPe - k14*CAR1)*dt
            
    
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

class Population:
    def __init__(self,pos,initialVals,params,t):
        self.pos=pos
        self.AgentParam=params
        self.t = t
        
        # preallocate variable arrays 
        self.ACA = np.zeros(len(t))
        self.PKA = np.zeros(len(t))
        self.ERK2 = np.zeros(len(t))
        self.RegA = np.zeros(len(t))
        self.cAMPi = np.zeros(len(t))
        self.cAMPe = np.zeros(len(t))
        self.CAR1 = np.zeros(len(t))
        
        # set initial values
        self.ACA[0] = initialVals[0]
        self.PKA[0] = initialVals[1]
        self.ERK2[0] = initialVals[2]
        self.RegA[0] = initialVals[3]
        self.cAMPi[0] = initialVals[4]
        self.cAMPe[0] = initialVals[5]
        self.CAR1[0] = initialVals[6] 
    
    def run(self,dt,cAMPe_in,rho=0,gamma=0):
        
        # get parameter values from dictionary
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
        
        # run simulation 
        for i in range(1,len(self.t)):
        
            # get variable values from previous timestep 
            ACA = self.ACA[i-1]
            PKA = self.PKA[i-1]
            ERK2 = self.ERK2[i-1]
            RegA = self.RegA[i-1]
            cAMPi = self.cAMPi[i-1]
            cAMPe = self.cAMPe[i-1]
            CAR1 = self.CAR1[i-1]
            
            # compute next values
            self.ACA[i] = ACA + (k1*CAR1 - k2*ACA*PKA)*dt
            self.PKA[i] = PKA + (k3*cAMPi - k4*PKA)*dt
            self.ERK2[i] = ERK2 + (k5*CAR1-k6*PKA*ERK2)*dt
            self.RegA[i] = RegA + (k7-k8*ERK2*RegA)*dt
            self.cAMPi[i] = cAMPi + (k9*ACA-k10*RegA*cAMPi)*dt
            self.cAMPe[i] = cAMPe + (k11*rho*ACA-(k12+gamma)*(cAMPe - cAMPe_in[i]))*dt
            self.CAR1[i] = CAR1 + (k13*cAMPe - k14*CAR1)*dt
            
            # correct cAMPe if less than 0
            if self.cAMPe[i] < 0 :
                self.cAMPe[i] = 0
           
        
    
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
        CAR1_next=self.CAR1_now+(k13*self.cAMPe_now- k14*self.CAR1_now)*dt
        
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

class Population_scNoise:
    def __init__(self,pos,initialVals,params,t):
        self.pos=pos
        self.AgentParam=params
        self.t = t
        N = self.AgentParam['N']
        
        # preallocate variable arrays 
        self.ACA = np.zeros([N,len(t)])
        self.PKA = np.zeros([N,len(t)])
        self.ERK2 = np.zeros([N,len(t)])
        self.RegA = np.zeros([N,len(t)])
        self.cAMPi = np.zeros([N,len(t)])
        self.cAMPe = np.zeros([len(t)])
        self.CAR1 = np.zeros([N,len(t)])
        
        # set initial values
        self.ACA[:,0] = initialVals[0]
        self.PKA[:,0] = initialVals[1]
        self.ERK2[:,0] = initialVals[2]
        self.RegA[:,0] = initialVals[3]
        self.cAMPi[:,0] = initialVals[4]
        self.cAMPe[0] = initialVals[5]
        self.CAR1[:,0] = initialVals[6] 
    
    def run(self,dt,cAMPe_in,r,rho=1,gamma=0):
        
        # get parameter values from dictionary
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
        N = self.AgentParam['N']
        
        # run simulation 
        for i in range(1,len(self.t)):
        
            # get variable values from previous timestep 
            ACA = self.ACA[:,i-1]
            PKA = self.PKA[:,i-1]
            ERK2 = self.ERK2[:,i-1]
            RegA = self.RegA[:,i-1]
            cAMPi = self.cAMPi[:,i-1]
            cAMPe = self.cAMPe[i-1]
            CAR1 = self.CAR1[:,i-1]
            
            # compute next values
            self.ACA[:,i] = ACA + (k1*CAR1 - k2*ACA*PKA)*dt
            self.PKA[:,i] = PKA + (k3*cAMPi - k4*PKA)*dt
            self.ERK2[:,i] = ERK2 + (k5*CAR1-k6*PKA*ERK2)*dt
            self.RegA[:,i] = RegA + (k7-k8*ERK2*RegA)*dt
            self.cAMPi[:,i] = cAMPi + (k9*ACA-k10*RegA*cAMPi)*dt + r[i]
            self.cAMPe[i] = cAMPe + (k11*rho*np.mean(ACA)-(k12+gamma)*(cAMPe - cAMPe_in[i]))*dt
            self.CAR1[:,i] = CAR1 + (k13*cAMPe - k14*CAR1)*dt
            
            # correct cAMPe if less than 0
            if self.cAMPe[i] < 0 :
                self.cAMPe[i] = 0
    
    def update(self,dt,campExt_influx,rho = 1,gamma = 0, r=np.array([])):
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
        
        N=self.AgentParam['N'] # number of cells in the population
        sigma = self.AgentParam['sigma'] # noise strength
        # don't fix random seed when r is empty array, fix random seed if r isn't a random array
        if r.size == 0: 
            r=  math.sqrt(dt)*np.random.normal(0,1,N)
        
        ACA_next=self.ACA_now+(k1*self.CAR1_now-k2*np.multiply(self.ACA_now, self.PKA_now))*dt
        PKA_next=self.PKA_now+(k3*self.cAMPi_now-k4*self.PKA_now)*dt
        ERK2_next=self.ERK2_now+(k5*self.CAR1_now-k6*np.multiply(self.PKA_now, self.ERK2_now))*dt
        RegA_next=self.RegA_now+(k7-k8*np.multiply(self.ERK2_now, self.RegA_now))*dt
        cAMPi_next=self.cAMPi_now+(k9*self.ACA_now-k10*np.multiply(self.RegA_now, self.cAMPi_now))*dt + sigma*r
        
        cAMPe_next=self.cAMPe_now+(k11*rho*np.sum(self.ACA_now)/N-(k12+gamma)*(self.cAMPe_now - campExt_influx))*dt
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