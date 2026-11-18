# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 19:23:30 2018

@author: Chuqiao
"""
import numpy as np
from numba import njit
from popHistory import CellPopulation, History, noise_generator, np_heaviside, put, stimulus

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
    
    
    def run(self,dt,time_separation,extracellularA_trace,rng=None,keep_cells=False):
        # rng is the random generator that the noise of each time step is drawn from. keep_cells
        # says for how many cells A and R are kept: every one (True), the first so many (a
        # number) or none (False). A_mean and R_mean are always kept. See popHistory.py

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
        
        # set initial A and R for all cells
        A = self.initialVals[0]
        R = self.initialVals[1]

        # what is kept of the run, and the stimulus and the noise as the compiled loop
        # takes them
        history = History(self,['A','R'],[A,R],keep_cells)
        extracellularA_trace = stimulus(extracellularA_trace,len(self.t))
        rng = noise_generator(rng,N)

        # take the time steps; the loop is steps(), below
        self.steps(A,R,e,g,c0,sigma,N,a,alpha0,alpha_pde,Kd,S,rho,j,time_separation,
                   self.cAMPext,extracellularA_trace,dt,self.t,rng,history.arrays)

        history.finish()

    @staticmethod
    @njit(error_model='numpy',cache=True)
    def steps(A,R,e,g,c0,sigma,N,a,alpha0,alpha_pde,Kd,S,rho,j,time_separation,self_cAMPext,
              extracellularA_trace,dt,t,rng,history):
        # The loop of run(). numba compiles this method: popHistory.py says what to keep
        # to when changing it.

        for i in range(1,len(t)):
            # draw the noise of this time step
            r = np.sqrt(dt)*rng.normal(0,1,N)

            # get cAMPe value from previous timestep; A and R are those of that timestep
            cAMPe = self_cAMPext[i-1]
            
            # compute differential 
            dA=A-(A**3.0)/3-R+a*np.log(1+cAMPe/Kd)
            dR=e*(A-g*R+c0)
        
            # calculate next external cAMP step with/without time separation
            alphaf = extracellularA_trace[i]
            if time_separation == 0:
                fcAMPext = alphaf + rho*alpha0 + rho*S/N*np.sum(np_heaviside(A,0.5))-(j+alpha_pde*rho) * cAMPe
                self_cAMPext[i] = cAMPe + fcAMPext*dt
            else:
                self_cAMPext[i] = (alphaf + rho*alpha0 + rho*S/N*np.sum(np_heaviside(A,0.5)))/(j+alpha_pde*rho)
            
            # compute time step 
            A = A + dA*dt + sigma*r
            R = R + dR*dt
            put(history,i,A,R)
