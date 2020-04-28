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
        
    
    def update(self,dt,signal):
        e=self.AgentParam['e']   
        g=self.AgentParam['g'] 
        c0=self.AgentParam['c0'] 
        sigma=self.AgentParam['sigma']
        # N=self.AgentParam['N'] 
        a=self.AgentParam['a'] 
        # alpha0=self.AgentParam['alpha0'] 
        # alpha_pde=self.AgentParam['alpha_pde'] 
        Kd=self.AgentParam['Kd'] 
        # S=self.AgentParam['S'] 
        
        
        fA=self.A_now-(self.A_now)**3/3-self.R_now+a*np.log(1+signal/Kd)
        fR=e*(self.A_now-g*self.R_now+c0)
        r= math.sqrt(dt)*np.random.normal(0, 1)
        A_next=self.A_now+fA*dt+sigma*r
        R_next=self.R_now+fR*dt
        
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
        
        
        
        
class Sgro2015_pop:
    def __init__(self,A0,R0,cAMPext0, PopParam):
        self.PopParam=PopParam
        self.A_now=A0
        self.R_now=R0 # initial state input as a list variable [A0,R0]
        self.cAMPext_now = cAMPext0
        
    
    def update(self, dt, time_separation, alphaf): # alphaf is exernal cAMP influx
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
        
        # r = math.sqrt(dt)*np.random.uniform(-1,1,N)
        r = math.sqrt(dt)*np.random.normal(0, 1, N)
        fA=self.A_now-(np.power(self.A_now,3))/3-self.R_now+a*np.log(1+self.cAMPext_now/Kd)
        fR=e*(self.A_now-g*self.R_now+c0)
        if time_separation == 0:
            fcAMPext = alphaf + rho*alpha0 + rho*S/N*np.sum(np.heaviside(self.A_now,0.5))-(j+alpha_pde*rho) * self.cAMPext_now
            cAMPext_next = self.cAMPext_now +dt*fcAMPext
        else:
            cAMPext_next = (alphaf + rho*alpha0 + rho*S/N*np.sum(np.heaviside(self.A_now,0.5)))/(j+alpha_pde*rho)
        
        A_next=self.A_now+fA*dt+sigma*r
        R_next=self.R_now+fR*dt
        
        
        self.A_now=A_next
        self.R_now=R_next # update the current A and R state
        self.cAMPext_now = cAMPext_next
        return A_next,R_next, cAMPext_next


    
    def flux(self,signals):
        if self.A>self.AgentParam['flux_thrs']:
            agent_flux=True
        else:
            agent_flux=False
        return agent_flux


    def print_state(self):
        print('past A:'+str(self.A_now))
        print('past R:'+str(self.R_now)) 


# choice ofcell-media coupling options 
class Sgro2015_pop_dir_cpl:
    def __init__(self,A0,R0,cAMPext0, PopParam):
        self.PopParam=PopParam
        self.A_now=A0
        self.R_now=R0 # initial state input as a list variable [A0,R0]
        self.cAMPext_now = cAMPext0
        
    
    def update(self, dt, time_separation, alphaf, dir_cpl): 
        '''
        alphaf: exernal cAMP influx (background cAMP)
        dir_cpl: whether and how are the cells directly coupled through external media
            dir_cpl = 0: no direct coupling, setup in the original model
            dir_cpl = 1: direct coupling through A
            dir_cpl = 2: direct coupling through A, with all A<0--> 0
            dir_cpl = 3: direct coupling through shifted, normalized A ((A+A_offset)/Nh_Sgro)
        '''
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
        
        # r = math.sqrt(dt)*np.random.uniform(-1,1,N)
        r = math.sqrt(dt)*np.random.normal(0, 1, N)
        fA=self.A_now-(np.power(self.A_now,3))/3-self.R_now+a*np.log(1+self.cAMPext_now/Kd)
        fR=e*(self.A_now-g*self .R_now+c0)
        
        if dir_cpl == 0:
            A_release = np.heaviside(self.A_now,0.5)
            
        elif dir_cpl == 1:
            temp = (self.A_now)*1; # direct coupling of raw value of A
            A_release = temp   
            
        elif dir_cpl == 2:
            temp = (self.A_now)*1; temp[temp<0] = 0 # prevent temp act as pointer of self.A_now
            A_release = temp # dir_cpl2, set all negative A to 0
            
        elif dir_cpl == 3:
            A_offset=1.5; Nh_Sgro = 3.5; temp = (self.A_now)*1; temp[temp<0] = 0
            A_release = (temp + A_offset)/Nh_Sgro # dir_cpl3
            
        if time_separation == 0:
            fcAMPext = alphaf + rho*alpha0 + rho*S/N*np.sum(A_release)-(j+alpha_pde*rho) * self.cAMPext_now
            cAMPext_next = self.cAMPext_now +dt*fcAMPext
        else:
            cAMPext_next = (alphaf + rho*alpha0 + rho*S/N*np.sum(A_release))/(j+alpha_pde*rho)
        
        A_next=self.A_now+fA*dt+sigma*r
        R_next=self.R_now+fR*dt  
    
        self.A_now=A_next
        self.R_now=R_next # update the current A and R state
        self.cAMPext_now = cAMPext_next
        
        return A_next,R_next, cAMPext_next


    
    def flux(self,signals):
        if self.A>self.AgentParam['flux_thrs']:
            agent_flux=True
        else:
            agent_flux=False
        return agent_flux


    def print_state(self):
        print('past A:'+str(self.A_now))
        print('past R:'+str(self.R_now)) 

class Sgro2015_pop_mixed_cells:
    def __init__(self,A0,R0,cAMPext0, PopParam):
        self.PopParam=PopParam
        self.A_now=A0
        self.R_now=R0 # initial state input as a list variable [A0,R0]
        self.cAMPext_now = cAMPext0
        
    
    def update(self, dt, time_separation, alphaf): # alphaf is exernal cAMP influx
        e=self.PopParam['e']  # e is a numpy array, can be set diffeently for each cell 
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
        
        # r = math.sqrt(dt)*np.random.uniform(-1,1,N)
        r = math.sqrt(dt)*np.random.normal(0, 1, N)
        fA=self.A_now-(np.power(self.A_now,3))/3-self.R_now+a*np.log(1+self.cAMPext_now/Kd)
        fR=np.multiply(e,(self.A_now-g*self.R_now+c0))
        if time_separation == 0:
            fcAMPext = alphaf + rho*alpha0 + rho*S/N*np.sum(np.heaviside(self.A_now,0.5))-(j+alpha_pde*rho) * self.cAMPext_now
            cAMPext_next = self.cAMPext_now +dt*fcAMPext
        else:
            cAMPext_next = (alphaf + rho*alpha0 + rho*S/N*np.sum(np.heaviside(self.A_now,0.5)))/(j+alpha_pde*rho)
        
        A_next=self.A_now+fA*dt+sigma*r
        R_next=self.R_now+fR*dt
        
        
        self.A_now=A_next
        self.R_now=R_next # update the current A and R state
        self.cAMPext_now = cAMPext_next
        return A_next,R_next, cAMPext_next

    def flux(self,signals):
        if self.A>self.AgentParam['flux_thrs']:
            agent_flux=True
        else:
            agent_flux=False
        return agent_flux
    def print_state(self):
        print('past A:'+str(self.A_now))
        print('past R:'+str(self.R_now))      
        
