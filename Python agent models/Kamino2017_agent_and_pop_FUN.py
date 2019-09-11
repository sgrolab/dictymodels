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

class Kamino2017_agent:
    def __init__(self,state,Param):
        self.state=state
        self.Param=Param
        self.x_now=state[0]
        self.y_now=state[1] 
        self.z_now=state[2] # initial state input as a list variable [x0,y0,z0]
        
    
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
class Kamino2017_pop:
    def __init__(self,state,Param):
        self.state=state
        self.Param=Param
        self.x_now=state[0]
        self.y_now=state[1] 
        self.z_now=state[2] # initial state input as a list variable [x0,y0,z0]
        
    
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

