# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 19:54:45 2019

@author: Chuqiao
"""
import numpy as np
from numba import njit
from popHistory import CellPopulation, History, noise_generator, np_sum, put, stimulus

class Cell():
    def __init__(self,pos,state,AgentParam,t):
        self.pos=pos
        self.AgentParam=AgentParam
        self.t = t
        
        # initialize variable arrays 
        self.cAMPi = np.zeros(len(t))
        self.thetai = np.zeros(len(t))
        self.cAMPe = np.zeros(len(t))
        
        # set initial values 
        self.cAMPi[0] = state[0]
        self.thetai[0] = state[1]
        self.cAMPe[0] = state[2]
        
    def run(self,dt,r,cAMPe_in):
        
        # pull parameters 
        Amax=self.AgentParam['Amax']   
        Abas=self.AgentParam['Abas']   
        w=self.AgentParam['w']   # min-1
        K=self.AgentParam['K']   
        c_excite=self.AgentParam['c_excite']
        
        # run simulation 
        for i in range(1,len(self.t)):
            
            # get values from pervious timestep
            thetai = self.thetai[i-1]
            
            # compute change in theta
            dthetai = w*(1-K/(cAMPe_in[i]+K)*c_excite*np.sin(thetai))
            
            # compute next time step 
            self.thetai[i] = thetai + dt*dthetai + r[i]
            self.cAMPi[i] = ((-Amax+Abas)*np.sin(self.thetai[i])+Amax+Abas)/2


class Gregor2010_pop(CellPopulation):
    def __init__(self,campCyto,thetai,campExt,PopParam,t):
        # get parameters 
        self.PopParam=PopParam
        self.t = t
        self.initialVals = [campCyto,thetai]
        
        # preallocate external cAMP matrix, which the cells share, and set its initial value
        self.cAMPe = np.zeros(len(t))
        self.cAMPe[0] = campExt
        
    def run(self,dt,eta,rho,k,Vt,time_separation,campExt_influx,rng=None,
            keep_cells=False):
        # rng is the random generator that the noise of each time step is drawn from. keep_cells
        # says for how many cells cAMPi and thetai are kept: every one (True), the first so many
        # (a number) or none (False). cAMPi_mean and thetai_mean are always kept. See
        # popHistory.py

        # pull model parameters
        Amax=self.PopParam['Amax']   
        Abas=self.PopParam['Abas']   
        w=self.PopParam['w']   # min-1
        Vc=self.PopParam['Vc']  
        St=self.PopParam['St']   
        Sc=self.PopParam['Sc']   
        K=self.PopParam['K']   
        c_sec=self.PopParam['c_sec']  
        c_excite=self.PopParam['c_excite']   
        Nc=self.PopParam['Nc']
        
        # set initial values, one for each cell
        cAMPi = np.zeros(Nc)
        thetai = np.zeros(Nc)
        cAMPi[:] = self.initialVals[0]
        thetai[:] = self.initialVals[1]

        # what is kept of the run, and the stimulus and the noise as the compiled loop
        # takes them
        history = History(self,['cAMPi','thetai'],[cAMPi,thetai],keep_cells)
        campExt_influx = stimulus(campExt_influx,len(self.t))
        rng = noise_generator(rng,Nc)

        # take the time steps; the loop is steps(), below
        self.steps(cAMPi,thetai,Amax,Abas,w,Vc,St,Sc,K,c_sec,c_excite,Nc,eta,rho,k,Vt,
                   time_separation,self.cAMPe,campExt_influx,dt,self.t,rng,history.arrays)

        history.finish()

    @staticmethod
    @njit(error_model='numpy',cache=True)
    def steps(cAMPi,thetai,Amax,Abas,w,Vc,St,Sc,K,c_sec,c_excite,Nc,eta,rho,k,Vt,
              time_separation,self_cAMPe,campExt_influx,dt,t,rng,history):
        # The loop of run(). numba compiles this method: popHistory.py says what to keep
        # to when changing it.
        
        # run model for full time 
        for i in range(1,len(t)):
            # draw the noise of this time step
            r = np.sqrt(dt)*rng.normal(0,1,Nc)

            # get cAMPe values from previous timestep
            cAMPe = self_cAMPe[i-1]
            cAMPe_in = campExt_influx[i]
    
            # calculate next external cAMP step with/without time separation
            if time_separation == 1:
                self_cAMPe[i] = Vc*St/Sc*rho/k*c_sec*1/Nc*np_sum(cAMPi) + cAMPe_in
            else:
                self_cAMPe[i] = cAMPe + dt*(Vc*St/Sc*rho/Vt*c_sec*1/Nc*np_sum(cAMPi)-k/Vt*(cAMPe-cAMPe_in))

            # calculate next time step for intracellular variables, both from thetai of
            # the previous timestep
            thetai, cAMPi = (
                thetai + dt*(w*(1-K/(cAMPe+K)*c_excite*np.sin(thetai))) + eta*r,
                ((-Amax+Abas)*np.sin(thetai)+(Amax+Abas))/2)
            put(history,i,cAMPi,thetai)
