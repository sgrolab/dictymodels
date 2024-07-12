# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 10:35:46 2019

@author: Chuqiao
"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import find_peaks

# remove all variables in the work space
# sys.modules[__name__].__dict__.clear()

class Cell:
    def __init__(self,initialVals,Param,t):
        self.t = t
        self.Param=Param
        
        # preallocate variable arrays 
        self.x = np.zeros(len(t))
        self.y = np.zeros(len(t))
        self.z = np.zeros(len(t))
        
        # set initial values 
        self.x[0] = initialVals[0]
        self.y[0] = initialVals[1]
        self.z[0] = initialVals[2]
    
    def run(self,dt,cAMPe_in):
        
        # pull parameter values 
        tau=self.Param['tau']   
        n=self.Param['n'] 
        K=self.Param['K'] 
        delta=self.Param['delta'] 
        
        # run simulation 
        for i in range(1,len(self.t)):
            
            # get variable values from previous timestep 
            x = self.x[i-1]
            y = self.y[i-1]
            z = self.z[i-1]
            
            # compute variable change 
            dx = (z+delta-x)/tau
            dy = (z+delta)**n/((z+delta)**n+(K*x)**n)-y
            
            # compute next variable values 
            self.x[i] = x + dx*dt
            self.y[i] = y + dy*dt
            self.z[i] = cAMPe_in[i]
            
    
    def update(self, dt, signal):
        tau=self.Param['tau']   
        n=self.Param['n'] 
        K=self.Param['K'] 
        delta=self.Param['delta'] 
#        kt=self.Param['kt'] 
#        gamma=self.Param['gamma'] 
#        rho=self.Param['rho'] 
        
        dxdt=((self.z_now+delta)-self.x_now)/tau
        dydt=(self.z_now+delta)**n/((self.z_now+delta)**n+(K*self.x_now)**n)-self.y_now
        
        x_next=self.x_now+dxdt*dt
        y_next=self.y_now+dydt*dt
        z_next = signal
        
        self.x_now=x_next
        self.y_now=y_next 
        self.z_now=z_next # update the current x,y,z state
        
        return x_next,y_next,z_next
    


    def print_state(self):
        print('past x:'+str(self.x_now))
        print('past y:'+str(self.y_now)) 
        print('past z:'+str(self.z_now)) 

# This is the same model with Kamino2017_agent      
class Population:
    def __init__(self,initialVals,params,t):
        self.Param=params
        self.t = t
        
        # preallocate variable arrays 
        self.x = np.zeros(len(t))
        self.y = np.zeros(len(t))
        self.z = np.zeros(len(t))
        
        # set initial values 
        self.x[0] = initialVals[0]
        self.y[0] = initialVals[1]
        self.z[0] = initialVals[2]
        
    def run(self,z0_influx, dt, overriding_sig='none'):
        
        # get parameter values
        tau=self.Param['tau']   
        n=self.Param['n'] 
        K=self.Param['K'] 
        delta=self.Param['delta'] 
        kt=self.Param['kt'] 
        gamma=self.Param['gamma'] 
        rho=self.Param['rho'] 
        
        # run simulation 
        for i in range(1,len(self.t)):
            
            # parameter values from last time step 
            x = self.x[i-1]
            y = self.y[i-1]
            z = self.z[i-1]
            
            # calculate next values 
            dx = (z + delta - x)/tau
            dy = (z + delta)**n/((z + delta)**n + (K*x)**n) - y
            dz = rho*kt*y - gamma*(z-z0_influx[i])
        
            # apply time step
            self.x[i] = x + dx*dt
            self.y[i] = y + dy*dt
            self.z[i] = z + dz*dt
        
    def update(self, z0_influx, dt, overriding_sig = 'none'):
        tau=self.Param['tau']   
        n=self.Param['n'] 
        K=self.Param['K'] 
        delta=self.Param['delta'] 
        kt=self.Param['kt'] 
        gamma=self.Param['gamma'] 
        rho=self.Param['rho'] 
        
        dxdt=((self.z_now+delta)-self.x_now)/tau
        dydt=(self.z_now+delta)**n/((self.z_now+delta)**n+(K*self.x_now)**n)-self.y_now
        
        x_next=self.x_now+dxdt*dt
        y_next=self.y_now+dydt*dt
        dzdt=rho*kt*self.y_now-gamma*(self.z_now-z0_influx)
        z_next=self.z_now+dzdt*dt
#        if isinstance(overriding_sig, str):
#            dzdt=rho*kt*self.y_now-gamma*(self.z_now-z0_influx)
#            z_next=self.z_now+dzdt*dt
#        else:
#            z_next = overriding_sig
        
        self.x_now=x_next
        self.y_now=y_next 
        self.z_now=z_next # update the current x,y,z state
        
        return x_next,y_next,z_next
    


    def print_state(self):
        print('past x:'+str(self.x_now))
        print('past y:'+str(self.y_now)) 
        print('past z:'+str(self.z_now)) 

class Population_scNoise: # population with noisy single cells 
    def __init__(self,initialVals,params,t):
        self.Param=params
        self.t = t
        N = self.Param['N']
        
        # preallocate variable arrays 
        self.x = np.zeros([N,len(t)])
        self.y = np.zeros([N,len(t)])
        self.z = np.zeros([N,len(t)])
        
        # set initial values 
        self.x[:,0] = initialVals[0]
        self.y[:,0] = initialVals[1]
        self.z[:,0] = initialVals[2]
    
    def run(self,z0_influx, dt, r, overriding_sig='none'):
        
        # get parameter values
        tau=self.Param['tau']   
        n=self.Param['n'] 
        K=self.Param['K'] 
        delta=self.Param['delta'] 
        kt=self.Param['kt'] 
        gamma=self.Param['gamma'] 
        rho=self.Param['rho'] 
        
        # run simulation 
        for i in range(1,len(self.t)):
            
            # parameter values from last time step 
            x = self.x[:,i-1]
            y = self.y[:,i-1]
            z = self.z[:,i-1]
            
            # calculate next values 
            dx = (z + delta - x)/tau
            dy = (z + delta)**n/((z + delta)**n + (K*x)**n) - y
            dz = rho*kt*y - gamma*(z-z0_influx[i])
        
            # apply time step
            self.x[:,i] = x + dx*dt
            self.y[:,i] = y + dy*dt + r[i]
            self.z[:,i] = z + dz*dt
    
    def update(self, z0_influx, dt,r,overriding_sig=None): 
        tau=self.Param['tau']   
        n=self.Param['n'] 
        K=self.Param['K'] 
        kt=self.Param['kt'] 
        gamma=self.Param['gamma'] 
        delta=self.Param['delta'] 
        rho=self.Param['rho'] 
        sigma=self.Param['sigma'] # sigma is noise strength in single cells
        N = self.Param['N'] # number of cells in the population
        
        if r.size==0:
            r = math.sqrt(dt)*np.random.normal(0, 1, N)  # random number for noise
            
        dxdt=((self.z_now+delta)-self.x_now)/tau
        dydt=np.power((self.z_now+delta),n)/(np.power((self.z_now+delta),n)+np.power((K*self.x_now),n))-self.y_now
        
        x_next=self.x_now+dxdt*dt 
        y_next=self.y_now+dydt*dt + sigma*r # noise added to y
        if not overriding_sig:
            dzdt=rho*kt*np.sum(self.y_now)/N-gamma*(self.z_now-z0_influx)
            z_next=self.z_now+dzdt*dt
        else:
            z_next = overriding_sig
        
        self.x_now=x_next
        self.y_now=y_next 
        self.z_now=z_next # update the current x,y,z state
        
        return x_next,y_next,z_next
   
    def print_state(self):
        print('past x:'+str(self.x_now))
        print('past y:'+str(self.y_now)) 
        print('past z:'+str(self.z_now)) 
