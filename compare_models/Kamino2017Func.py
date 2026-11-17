# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 10:35:46 2019

@author: Chuqiao
"""
import numpy as np
from popHistory import CellPopulation, History, noise_rows

class Cell():
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
            
class Population():
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

class Population_scNoise(CellPopulation): # population with noisy single cells
    def __init__(self,initialVals,params,t):
        self.Param=params
        self.t = t
        self.initialVals = initialVals
        
    def run(self,z0_influx, dt, r, overriding_sig='none', keep_cells=False):
        # r is the noise, an array or a function; keep_cells says whether x, y and z are
        # kept for every cell or only as x_mean, y_mean and z_mean: see popHistory.py
        
        # get parameter values
        tau=self.Param['tau']   
        n=self.Param['n'] 
        K=self.Param['K'] 
        delta=self.Param['delta'] 
        kt=self.Param['kt'] 
        gamma=self.Param['gamma'] 
        rho=self.Param['rho'] 
        N = self.Param['N']

        # set initial values, one for each cell
        x = np.zeros(N)
        y = np.zeros(N)
        z = np.zeros(N)
        x[:] = self.initialVals[0]
        y[:] = self.initialVals[1]
        z[:] = self.initialVals[2]

        # what is kept of the run, and the noise of each time step after the first
        history = History(self,['x','y','z'],[x,y,z],keep_cells)
        noise = noise_rows(r,len(self.t))
        next(noise)
        
        # run simulation 
        for i in range(1,len(self.t)):
            
            # calculate next values from x, y and z of the last time step
            dx = (z + delta - x)/tau
            dy = (z + delta)**n/((z + delta)**n + (K*x)**n) - y
            dz = rho*kt*np.mean(y) - gamma*(z-z0_influx[i])
        
            # apply time step
            x = x + dx*dt
            y = y + dy*dt + next(noise)
            z = z + dz*dt
            history.put(i,x,y,z)

        history.finish()
