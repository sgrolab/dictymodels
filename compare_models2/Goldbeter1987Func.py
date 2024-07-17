# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 15:39:14 2019

@author: Chuqiao
"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt

# Define 4 variable class
class Goldbeter1987_agent_4var:
    def __init__(self,pos,state,param):
        self.pos=pos
        self.state=state
        self.param=param
        self.p_now=state[0]
        self.a_now=state[1]
        self.b_now=state[2] 
        self.g_now=state[3] # initial state input as a list variable [p0, b0, g0]
        
    
    def update(self,dt,signals):
        k1=self.param['k1']   
        k2=self.param['k2']  
        L1=self.param['L1']  
        L2=self.param['L2']  
        c=self.param['c']  
        lamda=self.param['lamda']  
        theta=self.param['theta']  
        e=self.param['e']  
        q=self.param['q'] 
        sig=self.param['sig']  
        v=self.param['v']  
        k=self.param['k']  
        ki=self.param['ki']  
        kt=self.param['kt']  
        kc=self.param['kc']  
        h=self.param['h']  
        
        f1 = (k1+k2*self.g_now)/(1+self.g_now)
        f2 = (k1*L1+k2*L2*c*self.g_now)/(1+c*self.g_now)
        # original model
        Y= self.p_now*self.g_now/(1+ self.g_now)
        Ysq = (Y)**2
        PI = self.a_now*(lamda*theta + e*Ysq)/(1 + theta*self.a_now + (1 + self.a_now)*e*Ysq)
        
        # p_next=0.8
        p_next = self.p_now+dt*(-f1*self.p_now+f2*(1-self.p_now))
        a_next = self.a_now+dt*(v-k*self.a_now-sig*PI)
        b_next = self.b_now+dt*(q*sig*PI-(ki+kt)*self.b_now)
        
        if signals == 0:
            g_next = self.g_now+dt*(kt/h*self.b_now-kc*self.g_now)
        else: 
            g_next = signals # if there is contant cAMP input to the system
        
        self.p_now = p_next
        self.a_now = a_next
        self.b_now = b_next
        self.g_now = g_next
        return p_next,a_next,b_next,g_next


# Define 3 variable class agent
class Cell():
    def __init__(self,params,initialVals,t):
        
        self.param = params 
        self.t = t
        
        # preallocate variable arrays
        self.p = np.zeros(len(t))
        self.b = np.zeros(len(t))
        self.g = np.zeros(len(t))
        
        # set initial values 
        self.p[0] = initialVals[0]
        self.b[0] = initialVals[1]
        self.g[0] = initialVals[2]
    
    def run(self,dt,a,signal):
        # pull parameter values 
        k1=self.param['k1']   
        k2=self.param['k2']  
        L1=self.param['L1']  
        L2=self.param['L2']  
        c=self.param['c']  
        lamda=self.param['lamda']  
        theta=self.param['theta']  
        e=self.param['e']  
        q=self.param['q'] 
        sig=self.param['sig']  
        v=self.param['v']  
        k=self.param['k']  
        ki=self.param['ki']  
        kt=self.param['kt']  
        kc=self.param['kc']  
        h=self.param['h'] 
        
        # run simulation 
        for i in range(1,len(self.t)):
            
            # get variable values from previous timestep 
            p = self.p[i-1]
            b = self.b[i-1]
            g = self.g[i-1]
            
            # calculate change 
            f1 = (k1+k2*g)/(1+g)
            f2 = (k1*L1+k2*L2*c*g)/(1+c*g)
            Ysq = (p*g/(1+g))**2
            PI = a*(lamda*theta + e*Ysq)/(1 + theta*a + (1 + a)*e*Ysq)
            
            # calculate next time step 
            self.p[i] = p + dt*(-f1*p + f2*(1-p))
            self.b[i] = b + dt*(q*sig*PI-(ki+kt)*b)
            self.g[i] = signal[i] # set g to constant cAMP
        
    def update(self,dt, a, signal):
        k1=self.param['k1']   
        k2=self.param['k2']  
        L1=self.param['L1']  
        L2=self.param['L2']  
        c=self.param['c']  
        lamda=self.param['lamda']  
        theta=self.param['theta']  
        e=self.param['e']  
        q=self.param['q'] 
        sig=self.param['sig']  
        v=self.param['v']  
        k=self.param['k']  
        ki=self.param['ki']  
        kt=self.param['kt']  
        kc=self.param['kc']  
        h=self.param['h']  
        
        f1 = (k1+k2*self.g_now)/(1+self.g_now)
        f2 = (k1*L1+k2*L2*c*self.g_now)/(1+c*self.g_now)
        Ysq = (self.p_now*self.g_now/(1+ self.g_now))**2
        PI = a*(lamda*theta + e*Ysq)/(1 + theta*a + (1 + a)*e*Ysq)
        
        # p_next=0.8
        p_next = self.p_now+dt*(-f1*self.p_now+f2*(1-self.p_now))
        b_next = self.b_now+dt*(q*sig*PI-(ki+kt)*self.b_now)
        g_next = signal # contant cAMP applied to single cell
        
        self.p_now = p_next
        self.b_now = b_next
        self.g_now = g_next
        
        return p_next, b_next, g_next
    
class Population():
    def __init__(self,initialVals,params,t):
        self.param = params
        self.t = t
        
        self.p = np.zeros(len(t))
        self.a = np.zeros(len(t))
        self.b = np.zeros(len(t))
        self.g = np.zeros(len(t))
        
        # set initial values
        self.p[0] = initialVals[0]
        self.a[0] = initialVals[1]
        self.b[0] = initialVals[2]
        self.g[0] = initialVals[3]
    
    def run(self,dt,a,cAMPe_in):
        # pull parameter values
        k1=self.param['k1']   
        k2=self.param['k2']  
        L1=self.param['L1']  
        L2=self.param['L2']  
        c=self.param['c']  
        lamda=self.param['lamda']  
        theta=self.param['theta']  
        e=self.param['e']  
        q=self.param['q'] 
        sig=self.param['sig']  
        v=self.param['v']  
        k=self.param['k']  
        ki=self.param['ki']  
        kt=self.param['kt']  
        kc=self.param['kc']  
        h=self.param['h']  
        
        # run model for full time 
        for i in range(1,len(self.t)):
        
            # get variable values from previous timestep
            p = self.p[i-1]
            b = self.b[i-1]
            g = self.g[-1]
            
            # calculate model equations
            f1 = (k1+k2*g)/(1+g)
            f2 = (k1*L1+k2*L2*c*g)/(1+c*g)
            Ysq = (p*g/(1+g))**2
            PI = a*(lamda*theta + e*Ysq)/(1 + theta*a + (1 + a)*e*Ysq)
        
            # update variable values
            self.p[i] = p + dt*(-f1*p+f2*(1-p))
            self.b[i] = b + dt*(q*sig*PI-(ki+kt)*b)
            self.g[i] = g + dt*(kt/h*b-kc*(g-cAMPe_in[i])) # campExt_influx: contant cAMP input to the system
    
    def update(self,dt, a, campExt_influx):
        k1=self.param['k1']   
        k2=self.param['k2']  
        L1=self.param['L1']  
        L2=self.param['L2']  
        c=self.param['c']  
        lamda=self.param['lamda']  
        theta=self.param['theta']  
        e=self.param['e']  
        q=self.param['q'] 
        sig=self.param['sig']  
        v=self.param['v']  
        k=self.param['k']  
        ki=self.param['ki']  
        kt=self.param['kt']  
        kc=self.param['kc']  
        h=self.param['h']  
        
        f1 = (k1+k2*self.g_now)/(1+self.g_now)
        f2 = (k1*L1+k2*L2*c*self.g_now)/(1+c*self.g_now)
        Ysq = (self.p_now*self.g_now/(1+ self.g_now))**2
        PI = a*(lamda*theta + e*Ysq)/(1 + theta*a + (1 + a)*e*Ysq)
        
        # p_next=0.8
        p_next = self.p_now+dt*(-f1*self.p_now+f2*(1-self.p_now))
        b_next = self.b_now+dt*(q*sig*PI-(ki+kt)*self.b_now)
        g_next = self.g_now+dt*(kt/h*self.b_now-kc*(self.g_now-campExt_influx)) # campExt_influx: contant cAMP input to the system
        
        self.p_now = p_next
        self.b_now = b_next
        self.g_now = g_next
        
        return p_next, b_next, g_next
    
class Population_scNoise():
    def __init__(self,initialVals,params,t):
        self.param = params
        self.t = t
        N = self.param['N']
        
        self.p = np.zeros([N,len(t)])
        self.a = np.zeros([N,len(t)])
        self.b = np.zeros([N,len(t)])
        self.g = np.zeros([N,len(t)])
        
        # set initial values
        self.p[:,0] = initialVals[0]
        self.a[:,0] = initialVals[1]
        self.b[:,0] = initialVals[2]
        self.g[:,0] = initialVals[3]
    
    def run(self,dt,r,a,cAMPe_in):
        # pull parameter values 
        k1=self.param['k1']   
        k2=self.param['k2']  
        L1=self.param['L1']  
        L2=self.param['L2']  
        c=self.param['c']  
        lamda=self.param['lamda']  
        theta=self.param['theta']  
        e=self.param['e']  
        q=self.param['q'] 
        sig=self.param['sig']  
        v=self.param['v']  
        k=self.param['k']  
        ki=self.param['ki']  
        kt=self.param['kt']  
        kc=self.param['kc']  
        h=self.param['h']  
    
        # run model for full time 
        for i in range(1,len(self.t)):
        
            # get variable values from previous timestep
            p = self.p[:,i-1]
            b = self.b[:,i-1]
            g = self.g[:,i-1]
            
            # calculate model equations
            f1 = (k1 + k2*g)/(1+g)
            f2 = (k1*L1 + k2*L2*c*g)/(1+c*g)
            Ysq = (p*g/(1+g))**2
            PI = a*(lamda*theta + e*Ysq)/(1 + theta*a + (1 + a)*e*Ysq)
        
            # update variable values
            self.p[:,i] = p + dt*(-f1*p + f2*(1-p))
            self.b[:,i] = b + dt*(q*sig*PI - (ki+kt)*b) + r[i]
            self.g[:,i] = g + dt*(kt/h*np.mean(b)-kc*(g-cAMPe_in[i])) # campExt_influx: contant cAMP input to the system
        
    def update(self,dt, a, campExt_influx,noise):
        k1=self.param['k1']   
        k2=self.param['k2']  
        L1=self.param['L1']  
        L2=self.param['L2']  
        c=self.param['c']  
        lamda=self.param['lamda']  
        theta=self.param['theta']  
        e=self.param['e']  
        q=self.param['q'] 
        sig=self.param['sig']  
        v=self.param['v']  
        k=self.param['k']  
        ki=self.param['ki']  
        kt=self.param['kt']  
        kc=self.param['kc']  
        h=self.param['h']  
        
        f1 = (k1+k2*self.g_now)/(1+self.g_now)
        f2 = (k1*L1+k2*L2*c*self.g_now)/(1+c*self.g_now)
        Ysq = (self.p_now*self.g_now/(1+ self.g_now))**2
        PI = a*(lamda*theta + e*Ysq)/(1 + theta*a + (1 + a)*e*Ysq)
        r = math.sqrt(dt)*noise
        # p_next=0.8
        p_next = self.p_now+dt*(-f1*self.p_now+f2*(1-self.p_now))
        b_next = self.b_now+dt*(q*sig*PI-(ki+kt)*self.b_now)+r
        g_next = self.g_now+dt*(kt/h*self.b_now-kc*(self.g_now-campExt_influx)) # campExt_influx: contant cAMP input to the system
        
        self.p_now = p_next
        self.b_now = b_next
        self.g_now = g_next
        
        return p_next, b_next, g_next

class Goldbeter1987_pop_3var_SCnoise:
	def __init__(self,pos,p0,b0,g0,param):
		self.pos=pos
		self.param=param
		self.p_now=p0
		self.b_now=b0 
		self.g_now=g0 
	def update(self,dt, a, campExt_influx, r =np.array([])):
		k1=self.param['k1']   
		k2=self.param['k2']  
		L1=self.param['L1']  
		L2=self.param['L2']  
		c=self.param['c']  
		lamda=self.param['lamda']  
		theta=self.param['theta']  
		e=self.param['e']  
		q=self.param['q'] 
		sig=self.param['sig']  
		v=self.param['v']  
		k=self.param['k']  
		ki=self.param['ki']  
		kt=self.param['kt']  
		kc=self.param['kc']  
		h=self.param['h']  
		
		sigma = self.param['sigma']# single cell noise strength
		N = self.param['N']# number of cells in the population
		f1 = (k1+k2*self.g_now)/(1+self.g_now)
		f2 = (k1*L1+k2*L2*c*self.g_now)/(1+c*self.g_now)
		Ysq = (self.p_now*self.g_now/(1+ self.g_now))**2
		PI = a*(lamda*theta + e*Ysq)/(1 + theta*a + (1 + a)*e*Ysq)
		if r.size ==0:
        		r=math.sqrt(dt)*np.random.normal(0,1,N)
        
		p_next = self.p_now+dt*(-f1*self.p_now+f2*(1-self.p_now))
		b_next = self.b_now+dt*(q*sig*PI-(ki+kt)*self.b_now) +r*sigma # noise added to intracellular cAMP term
		g_next = self.g_now+dt*(kt/h*np.sum(self.b_now)/N-kc*(self.g_now-campExt_influx)) # campExt_influx: contant cAMP input to the system

		self.p_now = p_next
		self.b_now = b_next
		self.g_now = g_next
		
		return p_next, b_next, g_next