# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 19:23:30 2018

@author: Chuqiao
"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import random

class Cell:
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
    def __init__(self,A0,R0,cAMPext0,PopParam,t):
        self.PopParam=PopParam
        N=self.PopParam['N']
        self.t = t
        
        # preallocate A and R matrices         
        self.A = np.zeros([N,len(t)])
        self.R = np.zeros([N,len(t)])
        
        # preallocate external cAMP array 
        self.cAMPext = np.zeros(len(t))
        
        # set initial A and R for all cells 
        self.A[:,0] = A0*np.ones(N)
        self.R[:,0] = R0*np.ones(N) 
        
        # set initial cAMPe amount 
        self.cAMPext[0] = cAMPext0
    
    
    def run(self,dt,time_separation,extracellularA_trace,r=0):
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
        if r.size == 0:
            r = math.sqrt(dt)*np.random.normal(0,1,size = (len(self.t),N))
        
        for i in range(1,len(self.t)):
            # get A, R, and cAMPe values from previous timestep
            A = self.A[:,i-1]
            R = self.R[:,i-1]
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
            self.A[:,i] = A + dA*dt + sigma*r[i,:]
            self.R[:,i] = R + dR*dt
        
        
    def update(self, dt, time_separation, alphaf, r=np.array([])): # alphaf is exernal cAMP influx
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
         # don't fix random seed when r is empty array, fix random seed if r isn't a random array
        if r.size == 0:
            r=  math.sqrt(dt)*np.random.normal(0,1,N)
        
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
        
