# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 10:10:49 2019

@author: Chuqiao
"""

import numpy as np

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
