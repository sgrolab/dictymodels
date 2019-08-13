# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 19:23:30 2018

@author: Chuqiao
"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt

class Sgro2015_agent:
    def __init__(self,pos,state,AgentParam):
        self.pos=pos
        self.state=state
        self.AgentParam=AgentParam
        self.A_now=state[0]
        self.R_now=state[1] # initial state input as a list variable [A0,R0]
        
    
    def update(self,signals,dt):
        e=self.AgentParam['e']   
        tauA=self.AgentParam['tauA'] 
        tauR=self.AgentParam['tauR']
        g=self.AgentParam['g'] 
        c0=self.AgentParam['c0'] 
        sigma=self.AgentParam['sigma']
        N=self.AgentParam['N'] 
        a=self.AgentParam['a'] 
        alpha0=self.AgentParam['alpha0'] 
        alpha_pde=self.AgentParam['alpha_pde'] 
        Kd=self.AgentParam['Kd'] 
        S=self.AgentParam['S'] 
        
        
        fA=self.A_now-(self.A_now)**3/3-self.R_now+a*np.log(1+signals/Kd)
        fR=e*(self.A_now-g*self.R_now+c0)
        r= math.sqrt(dt)*random.uniform(-1,1)
        A_next=self.A_now+fA*dt+sigma*r
        R_next=self.R_now+fR*dt/tauR
        
        self.A_now=A_next
        self.R_now=R_next # update the current A and R state
        return A_next,R_next,r
    
    def flux(self,signals):
        if self.A>self.AgentParam['flux_thrs']:
            agent_flux=True
        else:
            agent_flux=False
        return agent_flux


    def print_state(self):
        print('past A:'+str(self.A_now))
        print('past R:'+str(self.R_now)) 


        