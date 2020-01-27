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
class Goldbeter1987_agent_3var(Goldbeter1987_agent_4var):
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
    
class Goldbeter1987_pop_3var(Goldbeter1987_agent_3var):
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