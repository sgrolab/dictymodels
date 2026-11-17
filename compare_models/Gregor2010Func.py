# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 19:54:45 2019

@author: Chuqiao
"""
import numpy as np
from popHistory import CellPopulation, History, noise_rows

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
        
    def run(self,dt,eta,rho,k,Vt,time_separation,campExt_influx,r =np.array([]),
            keep_cells=False):
        # r is the noise, an array or a function; keep_cells says whether cAMPi and thetai
        # are kept for every cell or only as cAMPi_mean and thetai_mean: see popHistory.py

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
        
        # if no random seed specified, generate random steps 
        if not callable(r) and r.size == 0:
            r = np.sqrt(dt)*np.random.normal(0,1,size = (len(self.t),Nc))

        # set initial values, one for each cell
        cAMPi = np.zeros(Nc)
        thetai = np.zeros(Nc)
        cAMPi[:] = self.initialVals[0]
        thetai[:] = self.initialVals[1]

        # what is kept of the run, and the noise of each time step after the first
        history = History(self,['cAMPi','thetai'],[cAMPi,thetai],keep_cells)
        noise = noise_rows(r,len(self.t))
        next(noise)
        
        # run model for full time 
        for i in range(1,len(self.t)):
            # get cAMPe values from previous timestep
            cAMPe = self.cAMPe[i-1]
            cAMPe_in = campExt_influx[i]
    
            # calculate next external cAMP step with/without time separation
            if time_separation == 1:
                self.cAMPe[i] = Vc*St/Sc*rho/k*c_sec*1/Nc*np.sum(cAMPi) + cAMPe_in
            else:
                self.cAMPe[i] = cAMPe + dt*(Vc*St/Sc*rho/Vt*c_sec*1/Nc*np.sum(cAMPi)-k/Vt*(cAMPe-cAMPe_in))

            # calculate next time step for intracellular variables, both from thetai of
            # the previous timestep
            thetai, cAMPi = (
                thetai + dt*(w*(1-K/(cAMPe+K)*c_excite*np.sin(thetai))) + eta*next(noise),
                ((-Amax+Abas)*np.sin(thetai)+(Amax+Abas))/2)
            history.put(i,cAMPi,thetai)

        history.finish()
