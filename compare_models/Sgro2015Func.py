# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 19:23:30 2018

@author: Chuqiao
"""
import numpy as np
from popHistory import CellPopulation, History, noise_rows

class Cell():
    def __init__(self,pos,state,AgentParam,t):
        self.pos=pos
        self.state=state
        self.AgentParam=AgentParam
        self.t = t 
        
        # preallocate variable arrays 
        self.A = np.zeros(len(t))
        self.R = np.zeros(len(t))
        
        # set initial values 
        self.A[0] = state[0]
        self.R[0] = state[1]
        
    def run(self,dt,signal,r = False):
        e=self.AgentParam['e']   
        g=self.AgentParam['g'] 
        c0=self.AgentParam['c0'] 
        sigma=self.AgentParam['sigma']
        a=self.AgentParam['a'] 
        Kd=self.AgentParam['Kd'] 
        
        for i in range(1,len(self.t)):
            
            # get values from previous time step 
            A = self.A[i-1]
            R = self.R[i-1]
            
            # compute change in variables 
            dA = A - A**3/3 - R + a*np.log(1+signal[i]/Kd)
            dR = e*(A-g*R+c0)
            
            # compute next value 
            self.A[i] = A + dA*dt + sigma*r[i]
            self.R[i] = R + dR*dt
    
        
class Sgro2015_pop(CellPopulation):
    def __init__(self,A0,R0,cAMPext0,PopParam,t):
        self.PopParam=PopParam
        N=self.PopParam['N']
        self.t = t
        
        # initial A and R for all cells
        self.initialVals = [A0*np.ones(N),R0*np.ones(N)]
        
        # preallocate external cAMP array, which the cells share
        self.cAMPext = np.zeros(len(t))
        
        # set initial cAMPe amount 
        self.cAMPext[0] = cAMPext0
    
    
    def run(self,dt,time_separation,extracellularA_trace,r=0,keep_cells=False):
        # r is the noise, an array or a function; keep_cells says whether A and R are kept
        # for every cell or only as A_mean and R_mean: see popHistory.py

        # pull parameters from dictionary 
        e=self.PopParam['e']   
        g=self.PopParam['g'] 
        c0=self.PopParam['c0'] 
        sigma=self.PopParam['sigma']
        N=self.PopParam['N'] 
        a=self.PopParam['a'] 
        alpha0=self.PopParam['alpha0'] 
        alpha_pde=self.PopParam['alpha_pde'] 
        Kd=self.PopParam['Kd'] 
        S=self.PopParam['S'] 
        rho =  self.PopParam['rho']
        j = self.PopParam['j']
        
        # if no random seed specified, generate 
        if not callable(r) and r.size == 0:
            r = np.sqrt(dt)*np.random.normal(0,1,size = (len(self.t),N))
        
        # set initial A and R for all cells
        A = self.initialVals[0]
        R = self.initialVals[1]

        # what is kept of the run, and the noise of each time step after the first
        history = History(self,['A','R'],[A,R],keep_cells)
        noise = noise_rows(r,len(self.t))
        next(noise)

        for i in range(1,len(self.t)):
            # get cAMPe value from previous timestep; A and R are those of that timestep
            cAMPe = self.cAMPext[i-1]
            
            # compute differential 
            dA=A-(A**3)/3-R+a*np.log(1+cAMPe/Kd)
            dR=e*(A-g*R+c0)
        
            # calculate next external cAMP step with/without time separation
            alphaf = extracellularA_trace[i]
            if time_separation == 0:
                fcAMPext = alphaf + rho*alpha0 + rho*S/N*np.sum(np.heaviside(A,0.5))-(j+alpha_pde*rho) * cAMPe
                self.cAMPext[i] = cAMPe + fcAMPext*dt
            else:
                self.cAMPext[i] = (alphaf + rho*alpha0 + rho*S/N*np.sum(np.heaviside(A,0.5)))/(j+alpha_pde*rho)
            
            # compute time step 
            A = A + dA*dt + sigma*next(noise)
            R = R + dR*dt
            history.put(i,A,R)

        history.finish()
