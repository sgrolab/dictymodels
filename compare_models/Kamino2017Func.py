# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 10:35:46 2019

@author: Chuqiao
"""
import numpy as np
from numba import njit
from popHistory import CellPopulation, History, noise_generator, np_mean, put, stimulus

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
        
    def run(self,z0_influx, dt, rng, overriding_sig='none', keep_cells=False, sigma=1):
        # rng is the random generator that the noise of each time step is drawn from, with
        # strength sigma. keep_cells says for how many cells x, y and z are kept: every one
        # (True), the first so many (a number) or none (False). x_mean, y_mean and z_mean are
        # always kept. See popHistory.py
        
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

        # what is kept of the run, and the stimulus and the noise as the compiled loop
        # takes them
        history = History(self,['x','y','z'],[x,y,z],keep_cells)
        z0_influx = stimulus(z0_influx,len(self.t))
        rng = noise_generator(rng,N)

        # take the time steps; the loop is steps(), below
        self.steps(x,y,z,tau,n,K,delta,kt,gamma,rho,z0_influx,dt,N,sigma,self.t,rng,
                   history.arrays)

        history.finish()

    @staticmethod
    @njit(error_model='numpy',cache=True)
    def steps(x,y,z,tau,n,K,delta,kt,gamma,rho,z0_influx,dt,N,sigma,t,rng,history):
        # The loop of run(). numba compiles this method: popHistory.py says what to keep
        # to when changing it.
        
        # run simulation 
        for i in range(1,len(t)):
            # draw the noise of this time step
            r = np.sqrt(dt)*rng.normal(0,1,N)*sigma
            
            # calculate next values from x, y and z of the last time step
            dx = (z + delta - x)/tau
            dy = (z + delta)**n/((z + delta)**n + (K*x)**n) - y
            dz = rho*kt*np_mean(y) - gamma*(z-z0_influx[i])
        
            # apply time step
            x = x + dx*dt
            y = y + dy*dt + r
            z = z + dz*dt
            put(history,i,x,y,z)
