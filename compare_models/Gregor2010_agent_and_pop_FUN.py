# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 19:54:45 2019

@author: Chuqiao
"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt

# import sys
# remove all variables in the work space
# sys.modules[__name__].__dict__.clear() 

class Gregor2010_agent:
    def __init__(self,pos,state,AgentParam):
        self.pos=pos
        self.state=state
        self.AgentParam=AgentParam
        self.campCyto_now=state[0]
        self.thetai_now = state[1]
        self.campExt_now=state[2] # initial state input as a list variable [A0,R0]
        
    
    def update(self,dt,eta, campExt ):
        Amax=self.AgentParam['Amax']   
        Abas=self.AgentParam['Abas']   
        w=self.AgentParam['w']   # min-1
        # Vc=self.AgentParam['Vc']  
        # St=self.AgentParam['St']   
        # Sc=self.AgentParam['Sc']   
        K=self.AgentParam['K']   
        # c_sec=self.AgentParam['c_sec']  
        c_excite=self.AgentParam['c_excite']   
        
        
        r=  math.sqrt(dt)*random.uniform(-1,1)
        
        thetai_next=self.thetai_now + dt*(w*(1-K/(campExt+K)*c_excite*math.sin(self.thetai_now))) + eta*r
        
        self.thetai_now=thetai_next    
        
        campCyto_next=((-Amax+Abas)*math.sin(thetai_next)+(Amax+Abas))/2
        
        self.campCyto_now=campCyto_next
        
        return thetai_next, campCyto_next, r
    
    
    def flux(self,signals):
        if self.A>self.AgentParam['flux_thrs']:
            agent_flux=True
        else:
            agent_flux=False
        return agent_flux


    def print_state(self):
        pass


class Gregor2010_pop:
    def __init__(self,campCyto,thetai,campExt,PopParam):
        self.PopParam=PopParam
        self.campCyto_now=campCyto
        self.thetai_now = thetai
        self.campExt_now=campExt 
        
    
    def update(self,dt,eta,rho,k,Vt,time_separation,campExt_influx):
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
        
        
        
        r=  math.sqrt(dt)*np.random.normal(0,1,Nc)
        
        thetai_next=self.thetai_now + dt*(w*(1-K/(self.campExt_now+K)*c_excite*np.sin(self.thetai_now))) + eta*r
        
        #self.thetai_now=thetai_next   
        
        campCyto_next=((-Amax+Abas)*np.sin(self.thetai_now)+(Amax+Abas))/2
        
        if time_separation == 1:
        
            campExt_next = Vc*St/Sc*rho/k*c_sec*1/Nc*np.sum(self.campCyto_now) + campExt_influx
            
        else:
            
            campExt_next = self.campExt_now + dt*(Vc*St/Sc*rho/Vt*c_sec*1/Nc*np.sum(self.campCyto_now)-k/Vt*(self.campExt_now-campExt_influx))
        
        self.campCyto_now=campCyto_next
        self.campExt_now=campExt_next
        self.thetai_now=thetai_next                                     
        return thetai_next, campCyto_next, campExt_next
    
    
    def flux(self,signals):
        if self.A>self.AgentParam['flux_thrs']:
            agent_flux=True
        else:
            agent_flux=False
        return agent_flux


    def print_state(self):
        pass



