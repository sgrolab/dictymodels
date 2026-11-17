# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 10:10:49 2019

@author: Chuqiao
"""

import numpy as np
from numba import njit
from popHistory import CellPopulation, History, noise_generator, np_mean, put, stimulus

class Cell():
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
            
class Population():
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

class Population_scNoise(CellPopulation):
    def __init__(self,pos,initialVals,params,t):
        self.pos=pos
        self.AgentParam=params
        self.t = t
        self.initialVals = initialVals
        
        # preallocate the array of external cAMP, which the cells share, and set its
        # initial value
        self.cAMPe = np.zeros([len(t)])
        self.cAMPe[0] = initialVals[5]

    def run(self,dt,cAMPe_in,rng,rho=1,gamma=0,keep_cells=False,sigma=1):
        # rng is the random generator that the noise of each time step is drawn from, with
        # strength sigma. keep_cells says whether ACA, PKA, ERK2, RegA, cAMPi and CAR1 are
        # kept for every cell or only as ACA_mean and so on. See popHistory.py
        
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

        # set initial values, one for each cell
        ACA = np.zeros(N)
        PKA = np.zeros(N)
        ERK2 = np.zeros(N)
        RegA = np.zeros(N)
        cAMPi = np.zeros(N)
        CAR1 = np.zeros(N)
        ACA[:] = self.initialVals[0]
        PKA[:] = self.initialVals[1]
        ERK2[:] = self.initialVals[2]
        RegA[:] = self.initialVals[3]
        cAMPi[:] = self.initialVals[4]
        CAR1[:] = self.initialVals[6]

        # what is kept of the run, and the stimulus and the noise as the compiled loop
        # takes them
        history = History(self,['ACA','PKA','ERK2','RegA','cAMPi','CAR1'],
                          [ACA,PKA,ERK2,RegA,cAMPi,CAR1],keep_cells)
        cAMPe_in = stimulus(cAMPe_in,len(self.t))
        rng = noise_generator(rng,N)

        # take the time steps; the loop is steps(), below
        self.steps(ACA,PKA,ERK2,RegA,cAMPi,CAR1,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,
                   rho,gamma,self.cAMPe,cAMPe_in,dt,N,sigma,self.t,rng,history.arrays)

        history.finish()

    @staticmethod
    @njit(error_model='numpy',cache=True)
    def steps(ACA,PKA,ERK2,RegA,cAMPi,CAR1,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,rho,
              gamma,self_cAMPe,cAMPe_in,dt,N,sigma,t,rng,history):
        # The loop of run(). numba compiles this method: popHistory.py says what to keep
        # to when changing it.
        
        # run simulation 
        for i in range(1,len(t)):
            # draw the noise of this time step
            r = np.sqrt(dt)*rng.normal(0,1,N)*sigma
        
            # get external cAMP from previous timestep
            cAMPe = self_cAMPe[i-1]
            
            # compute next values, all from the values of the previous timestep
            self_cAMPe[i] = cAMPe + (k11*rho*np_mean(ACA)-(k12+gamma)*(cAMPe - cAMPe_in[i]))*dt
            ACA, PKA, ERK2, RegA, cAMPi, CAR1 = (
                ACA + (k1*CAR1 - k2*ACA*PKA)*dt,
                PKA + (k3*cAMPi - k4*PKA)*dt,
                ERK2 + (k5*CAR1-k6*PKA*ERK2)*dt,
                RegA + (k7-k8*ERK2*RegA)*dt,
                cAMPi + (k9*ACA-k10*RegA*cAMPi)*dt + r,
                CAR1 + (k13*cAMPe - k14*CAR1)*dt)
            put(history,i,ACA,PKA,ERK2,RegA,cAMPi,CAR1)
            
            # correct cAMPe if less than 0
            if self_cAMPe[i] < 0 :
                self_cAMPe[i] = 0
