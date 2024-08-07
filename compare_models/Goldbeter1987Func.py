# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 15:39:14 2019

@author: Chuqiao
"""
import numpy as np

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
        ki=self.param['ki']  
        kt=self.param['kt']  
        
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
        