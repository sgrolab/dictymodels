# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 19:54:45 2019

@author: Chuqiao
"""
import numpy as np

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


class Gregor2010_pop():
    def __init__(self,campCyto,thetai,campExt,PopParam,t):
        # get parameters 
        self.PopParam=PopParam
        N = self.PopParam['Nc']
        self.t = t
        
        # preallocate variable matrices
        self.cAMPi = np.zeros([N,len(t)])
        self.thetai = np.zeros([N,len(t)])
        
        # preallocate external cAMP matrix
        self.cAMPe = np.zeros(len(t))
        
        # set initial values
        self.cAMPi[:,0] = campCyto
        self.thetai[:,0] = thetai
        self.cAMPe[0] = campExt
        
    def run(self,dt,eta,rho,k,Vt,time_separation,campExt_influx,r =np.array([])):
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
        if r.size == 0:
            r = np.sqrt(dt)*np.random.normal(0,1,size = (len(self.t),Nc))
        
        # run model for full time 
        for i in range(1,len(self.t)):
            # get cAMPi, thetai, and cAMPe values from previous timestep
            thetai = self.thetai[:,i-1]
            cAMPi = self.cAMPi[:,i-1]
            cAMPe = self.cAMPe[i-1]
            cAMPe_in = campExt_influx[i]
            
            # calculate next time step for intracellular variables 
            self.thetai[:,i] = thetai + dt*(w*(1-K/(cAMPe+K)*c_excite*np.sin(thetai))) + eta*r[i,:]
            self.cAMPi[:,i] = ((-Amax+Abas)*np.sin(thetai)+(Amax+Abas))/2
    
            # calculate next external cAMP step with/without time separation
            if time_separation == 1:
                self.cAMPe[i] = Vc*St/Sc*rho/k*c_sec*1/Nc*np.sum(cAMPi) + cAMPe_in
            else:
                self.cAMPe[i] = cAMPe + dt*(Vc*St/Sc*rho/Vt*c_sec*1/Nc*np.sum(cAMPi)-k/Vt*(cAMPe-cAMPe_in))
