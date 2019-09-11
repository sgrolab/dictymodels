#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 17:01:34 2019

@author: goldalex
"""

from __future__ import division
import numpy as np
import scipy as sp
import time
import matplotlib.pyplot as plt
import scipy.sparse as sparse
import matplotlib.animation as animation
from randomgen import RandomGenerator
import line_profiler

profile = line_profiler.LineProfiler()

rnd = RandomGenerator()

cell_params_default={
        'c0' : 0.22,
        'a' : .058,
        'gamma' : 0.5,
        'Kd' : 10**(-5),
        'sigma' : .15,
        'epsilon' : 0.2,
        'cStim' : 100,
        'J' : 10,
        'rho' :  1,
        'D' : 1000,
        'a0' : 1,
        'af' : 0
       }

grid_params_default = {
        'dx' : 0.1,
        'D_sig' : 1,
        'box_size' : 50,
        'agent_dim' : .5,
        'p_vacant' : 0./100
        }

grid_params_default_masked = {
        'dx' : 0.1,
        'D_sig' : 1,
        'box_size_x' : 50,
        'box_size_y' : 50,
        'agent_dim' : .5,
        'num_agents' : 100*100,
        }


def accumulate_arr(coords, arr, shape):

    # Get linear indices to be used as IDs with bincount
    lidx = np.ravel_multi_index(coords, shape)
    # Or lidx = coords[0]*(coords[1].max()+1) + coords[1]

    # Accumulate arr with IDs from lidx
    return np.bincount(lidx,arr,minlength=shape[0]*shape[1]).reshape(shape[0],shape[1])

def scale(A, B, k):     # fill A with B scaled by k
    Y = A.shape[0]
    X = A.shape[1]
    for y in range(0, k):
        for x in range(0, k):
            A[y:Y:k, x:X:k] = B
    return A

def scale1D(A, B, k):     # fill A with B scaled by k
    Y = A.shape[0]
    for y in range(0, k):
        A[y:Y:k] = B
    return A


def scale_down_mat_sp(current_shape,scale):
    assert current_shape[0]//scale == current_shape[0]/scale, "Shape not divisible by scale"
    assert current_shape[1]//scale == current_shape[1]/scale, "Shape not divisible by scale"
    
    out_shape = (current_shape[0]//scale,current_shape[1]//scale)
    
    line = np.array(([1/scale**2]*scale+[0]*(current_shape[0]-scale))*scale+[0]*(current_shape[0]*(current_shape[1]-scale)))
    out_mat = sp.sparse.bsr_matrix(line)
    for i in range(out_shape[0]-1):
        newline = np.roll(line,scale*(i+1))
        out_mat = sp.sparse.vstack([out_mat,newline])
    block = out_mat.todense()
    for i in range(out_shape[1]-1):
        newblock = sp.sparse.bsr_matrix(np.roll(block,scale*current_shape[0]*(i+1),1))
        out_mat = sp.sparse.vstack([out_mat,newblock])
    return sparse.csr_matrix(out_mat)

def calc_square_laplacian_noflux_matrix(size):
    #set up basic laplacian block diagonal matrix
    xDiff = np.zeros((size[1]+1,size[1]))
    ix,jx = np.indices(xDiff.shape)
    xDiff[ix==jx] = 1
    xDiff[ix==jx+1] = -1
    
    yDiff = np.zeros((size[0]+1,size[0]))
    iy,jy = np.indices(yDiff.shape)
    yDiff[iy==jy] = 1
    yDiff[iy==jy+1] = -1
    
    Ax = sparse.dia_matrix(-np.matmul(np.transpose(xDiff),xDiff))
    Ay = sparse.dia_matrix(-np.matmul(np.transpose(yDiff),yDiff))
    
    lap =  sparse.kron(Ay,sparse.eye(size[1]),format='csr') + sparse.kron(sparse.eye(size[0]),Ax,format='csr')
    
    #set up boundary conditions
    
    lap += sparse.diags([2]+[1]*(size[1]-2)+[2]+([1]+[0]*(size[1]-2)+[1])*(size[0]-2)+[2]+[1]*(size[1]-2)+[2])
    lap = sparse.csr_matrix(lap)
    return lap

def lapMatGivenFlux(size,fluxes):
    
    
    return False

class dictyPop:
    
    def __init__(self,T,save_every,g_params=grid_params_default_masked,c_params=cell_params_default, noise = False,progress_report = True):
        self.g_params = g_params
        self.c_params = c_params
        self.T = T
        self.t = 0
        self.dt=self.g_params['dx']**2/(8*self.g_params['D_sig'])
        self.Tsteps=np.int(self.T/self.dt)
        self.signal_size = (self.g_params['box_size_x']/self.g_params['dx'],self.g_params['box_size_y']/self.g_params['dx']) #np.shape(np.eye(int(self.g_params['box_size']/self.g_params['dx'])))
        self.fluxes = np.zeros(4)
        self.p_r = progress_report
        self.save_every = save_every
        
        self.agent_ratio = int(self.g_params['agent_dim']/self.g_params['dx'])
        self.agent_size=(max(1,int(self.g_params['box_size_x']/self.g_params['agent_dim'])),max(1,int(self.g_params['box_size_y']/self.g_params['agent_dim'])))
        
        self.signal = np.zeros(self.signal_size)
        
        self.noise_flag = noise
        
        cell_state = rnd.normal(loc = 0.0, scale = 2, size = (2,g_params['num_agents']))
        
        self.A = cell_state[0,:]
        self.R = cell_state[1,:]
        
        self.dsignal = np.zeros(self.signal.shape)
        self.dA = np.zeros(self.A.shape)
        self.dR = np.zeros(self.R.shape)
        
        self.interp_mat = scale_down_mat_sp(self.signal_size,self.agent_ratio)
        self.A_grid = np.zeros((self.signal_size[0]//self.agent_ratio,self.signal_size[1]//self.agent_ratio))
        self.A_grid_big = np.zeros(self.signal.shape)
        
        self.coords = np.zeros((2,self.g_params['num_agents']),dtype = 'int')
        for agent in range(self.g_params['num_agents']):
            self.coords[:,agent] = np.array([agent%self.agent_size[0],agent//self.agent_size[0]])
            
        
        if min(self.signal_size) > 1:
            self.lap_mat = calc_square_laplacian_noflux_matrix(self.signal_size)
        else:
            self.lap_mat = sparse.diags([np.ones(max(self.signal_size)-1),-2*np.ones(max(self.signal_size)),np.ones(max(self.signal_size)-1)],[-1,0,1],format = 'csr')
            self.lap_mat += sparse.diags(np.array([1]+[0]*(max(self.signal_size)-2)+[1]),format = 'csr')
        
        self.A_saved = np.zeros((self.signal_size[0]//self.agent_ratio,self.signal_size[1]//self.agent_ratio,(save_every//self.T)+1))
        self.R_saved = np.zeros((self.signal_size[0]//self.agent_ratio,self.signal_size[1]//self.agent_ratio,(save_every//self.T)+1))
    
    def getAGrid(self):
        
        return accumulate_arr(self.coords,self.A,self.A_grid.shape)
    
    def getRGrid(self):
        
        return accumulate_arr(self.coords,self.R,self.A_grid.shape)
        
    def setSignal(self,signal):
        assert (self.signal.shape is signal.shape) and (self.signal.dtype is signal.dtype), "Input signal incorrect shape or type"
        self.signal = np.array(signal)
    
    def setCellState(self,state):
        assert (self.A.shape == state[:,0].shape) and (self.A.dtype == state.dtype), "Input cell state incorrect shape or type"
        self.A = np.array(state[:,0])
        self.R = np.array(state[:,1])
    
    def setCoords(self,coords):
        assert (self.coords.shape is coords[:,:,0].shape) and (self.coords.dtype is coords.dtype), "Input coordinates incorrect shape or type"
        self.coords = np.array(coords)
    
    def setFluxes(self,fluxes):
        assert (self.fluxes.shape is fluxes.shape) and (self.fluxes.dtype is fluxes.dtype), "Input fluxes incorrect shape or type"
        self.fluxes = np.array(fluxes)
        
    #@profile
    def getdA(self):
        return ((self.A-(self.A*self.A*self.A)/3 - self.R ))*self.dt\
        + ((self.c_params['a']*np.log1p(np.reshape(self.interp_mat.dot(np.reshape(self.signal,self.signal_size[0]*self.signal_size[1])),self.agent_size)[self.coords[0,:],self.coords[1,:]]/self.c_params['Kd'])))*self.dt 
        
    #@profile
    def getdR(self):
        
        return (((self.A-self.c_params['gamma']*self.R)+self.c_params['c0'])*self.c_params['epsilon'])*self.dt
    
    #@profile
    def getdC(self):
        
        #get laplacian for 0 flux
        laplacian = 1/(self.g_params['dx']**2)*np.reshape(self.lap_mat.dot(np.reshape(self.signal,self.signal_size[0]*self.signal_size[1])),(self.signal_size[0],self.signal_size[1]))
        
        #set fluxes
        
        self.signal[0,1:-1] -= self.fluxes[0]/self.g_params['dx']
        self.signal[-1,1:-1] -= self.fluxes[1]/self.g_params['dx']
        self.signal[1:-1,0] -= self.fluxes[2]/self.g_params['dx']
        self.signal[1:-1,-1] -= self.fluxes[3]/self.g_params['dx']
        
        self.signal[0,0] -= (self.fluxes[0]+self.fluxes[2])/(2*self.g_params['dx'])
        self.signal[-1,0] -= (self.fluxes[1]+self.fluxes[2])/(2*self.g_params['dx'])
        self.signal[0,-1] -= (self.fluxes[0]+self.fluxes[3])/(2*self.g_params['dx'])
        self.signal[-1,-1] -= (self.fluxes[1]+self.fluxes[3])/(2*self.g_params['dx'])
        
        self.A_grid  = accumulate_arr(self.coords,self.A,self.A_grid.shape)
        
        self.A_grid_big = scale(self.A_grid_big,self.A_grid,self.agent_ratio)
        
        return self.dt*(self.g_params['D_sig']*laplacian - self.c_params['J']*self.signal + self.c_params['rho']*self.c_params['a0'])\
        + self.dt*(self.c_params['D']*self.c_params['rho']*np.heaviside(self.A_grid_big,0.5))
        
    #@profile
    def update(self):
        
        self.dA = self.getdA()
        
        if self.noise_flag:
            self.dA += self.c_params['sigma']*rnd.normal(loc = 0.0, scale = np.sqrt(self.dt), size = self.A.shape)
        
        self.dR = self.getdR()
        
        self.dsignal = self.getdC()
        
        self.A += self.dA
        self.R += self.dR
        self.signal += self.dsignal
    
    #@profile
    def run(self):
        tic = time.time()
        j = 0
        for i in range(self.Tsteps):
            self.t += self.dt
            
            self.update()
            
            if ((i % (self.Tsteps//10)) == 0) and self.p_r == True:
                print((100*i/self.Tsteps),'% done')
                print('Time elapsed: ',time.time()-tic)
            
            if (i % (self.save_every//self.dt)) == 0:
                self.A_saved[:,:,j] = self.getAGrid()
                self.R_saved[:,:,j] = self.getRGrid()
                j += 1
                
            
        
class dictyPopOptoActive(dictyPop):
    def __init__(self,T,save_every,g_params=grid_params_default_masked,c_params=cell_params_default, noise = False,progress_report = True):
        self.g_params = g_params
        self.c_params = c_params
        self.T = T
        self.t = 0
        self.dt=self.g_params['dx']**2/(8*self.g_params['D_sig'])
        self.Tsteps=np.int(self.T/self.dt)
        self.signal_size=(int(self.g_params['box_size_x']/self.g_params['dx']),int(self.g_params['box_size_y']/self.g_params['dx']))
        self.fluxes = np.zeros(4)
        self.p_r = progress_report
        self.save_every = save_every
        
        self.agent_ratio = int(self.g_params['agent_dim']/self.g_params['dx'])
        self.agent_size=(max(1,int(self.g_params['box_size_x']/self.g_params['agent_dim'])),max(1,int(self.g_params['box_size_y']/self.g_params['agent_dim'])))
        
        self.signal = np.zeros(self.signal_size)
        
        self.noise_flag = noise
        
        cell_state = rnd.normal(loc = 0.0, scale = 2, size = (2,g_params['num_agents']))
        
        self.A = cell_state[0,:]
        self.R = cell_state[1,:]
        
        self.dsignal = np.zeros(self.signal.shape)
        self.dA = np.zeros(self.A.shape)
        self.dR = np.zeros(self.R.shape)
        
        self.interp_mat = scale_down_mat_sp(self.signal_size,self.agent_ratio)
        self.A_grid = np.zeros((self.signal_size[0]//self.agent_ratio,self.signal_size[1]//self.agent_ratio))
        self.A_grid_big = np.zeros(self.signal.shape)
        
        self.coords = np.zeros((2,self.g_params['num_agents']),dtype = 'int')
        for agent in range(self.g_params['num_agents']):
            self.coords[:,agent] = np.array([agent%self.agent_size[0],agent//self.agent_size[0]])
        
        self.mask = np.ones(self.g_params['num_agents'],dtype = 'int')
            

        
        self.lap_mat = calc_square_laplacian_noflux_matrix(self.signal_size)
    
        self.A_saved = np.zeros((self.signal_size[0]//self.agent_ratio,self.signal_size[1]//self.agent_ratio,(save_every//self.T)+1))
        self.R_saved = np.zeros((self.signal_size[0]//self.agent_ratio,self.signal_size[1]//self.agent_ratio,(save_every//self.T)+1))
    
    def update(self,mask_input_A,mask_input_R):
        
        self.dA = np.multiply(self.mask[np.newaxis,:],self.getdA())[0,:]
        
        if self.noise_flag:
            self.dA += np.multiply(self.mask[np.newaxis,:],self.c_params['sigma']*rnd.normal(loc = 0.0, scale = np.sqrt(self.dt), size = self.A.shape))[0,:]
        
        self.dR = np.multiply(self.mask[np.newaxis,:],self.getdR())[0,:]
        
        self.dsignal = self.getdC()
        
        self.A += self.dA
        self.R += self.dR
        self.signal += self.dsignal
        
        mask_coords = np.nonzero(1-self.mask)
        
        self.A[mask_coords] = mask_input_A
        
        self.R[mask_coords] = mask_input_R
        
        
    def run(self,mask_input_A,mask_input_R):
        tic = time.time()
        assert mask_input_A.size == self.Tsteps, 'input has incorrect length'
        assert mask_input_R.size == self.Tsteps, 'input has incorrect length'
        j = 0
        for i in range(self.Tsteps):
            self.t += self.dt
            self.update(mask_input_A[i],mask_input_R[i])
            
            if ((i % (self.Tsteps//10)) == 0) and self.p_r == True:
                print((100*i/self.Tsteps),'% done')
                print('Time elapsed: ',time.time()-tic)
            
            if (i % (self.save_every//self.dt)) == 0:
                self.A_saved[:,:,j] = self.getAGrid()
                self.R_saved[:,:,j] = self.getRGrid()
                j += 1
        
    
    def setMask(self,mask):
        assert (self.mask.shape == mask.shape) and (self.mask.dtype is mask.dtype), "Input mask incorrect shape or type"
        self.mask = mask
        
class dicty_pop_opto_active_1D(dictyPopOptoActive):
    def __init__(self,T,save_every,g_params=grid_params_default_masked,c_params=cell_params_default, noise = False,progress_report = True):
        self.g_params = g_params
        self.c_params = c_params
        self.T = T
        self.t = 0
        self.dt=self.g_params['dx']**2/(8*self.g_params['D_sig'])
        self.Tsteps=np.int(np.rint(self.T/self.dt))
        self.signal_size=int(self.g_params['box_size_x']/self.g_params['dx'])
        self.fluxes = np.zeros(2)
        self.p_r = progress_report
        self.save_every = save_every
        
        self.agent_ratio = int(self.g_params['agent_dim']/self.g_params['dx'])
        self.agent_size=int(self.g_params['box_size_x']/self.g_params['agent_dim'])
        
        self.signal = np.zeros(self.signal_size)
        
        self.noise_flag = noise
        
#        cell_state = -1*np.ones((g_params['num_agents'],2))
        
        self.A = -1*np.ones(g_params['num_agents'])
        self.R = np.zeros(g_params['num_agents'])
        
        self.dsignal = np.zeros(self.signal.shape)
        self.dA = np.zeros(self.A.shape)
        self.dR = np.zeros(self.R.shape)
        #get interpolation matrix
        line = np.zeros(self.signal_size)
        line[:self.agent_ratio] = 1/self.agent_ratio
        self.interp_mat = np.array(line)
        for i in range(self.agent_size-1):
            newline = np.roll(line,self.agent_ratio*(i+1))
            self.interp_mat = np.vstack((self.interp_mat,newline))
        
        self.A_grid = np.zeros((self.signal_size//self.agent_ratio))
        self.A_grid_big = np.zeros(self.signal.shape)
        
        self.coords = np.zeros(self.g_params['num_agents'],dtype = 'int')
        for agent in range(self.g_params['num_agents']):
            self.coords[agent] = agent
        
        self.mask = np.ones(self.g_params['num_agents'],dtype = 'int')
            

        
        self.lap_mat = sp.sparse.diags([np.ones(self.signal_size-1),-2*np.ones(self.signal_size),np.ones(self.signal_size-1)],[-1,0,1],format = 'csr')
        self.lap_mat += sp.sparse.diags(np.array([1]+[0]*(self.signal_size-2)+[1]),format = 'csr')
    
        self.A_saved = np.zeros((self.signal_size//self.agent_ratio,(self.T//save_every)+1))
        self.R_saved = np.zeros((self.signal_size//self.agent_ratio,(self.T//save_every)+1))
        
    def getdA(self):
        return ((self.A-(self.A*self.A*self.A)/3 - self.R ))*self.dt\
        + ((self.c_params['a']*np.log1p(self.interp_mat.dot(self.signal)[self.coords]/self.c_params['Kd'])))*self.dt
    
    #@profile
    def getdC(self):
        
        #get laplacian for 0 flux
        laplacian = 1/(self.g_params['dx']**2)*self.lap_mat.dot(self.signal)
        
        #set fluxes
        
        self.signal[0] -= self.fluxes[0]/self.g_params['dx']
        self.signal[-1] -= self.fluxes[1]/self.g_params['dx']
        
        self.A_grid  = accumulate_arr((np.zeros(self.coords.size,dtype = int),self.coords),self.A,(self.A_grid.shape[0],self.A_grid.shape[0]))[0,:]
        
        self.A_grid_big = scale1D(self.A_grid_big,self.A_grid,self.agent_ratio)
        
        return self.dt*(self.g_params['D_sig']*laplacian - self.c_params['J']*self.signal + self.c_params['rho']*self.c_params['a0'])\
        + self.dt*(self.c_params['D']*self.c_params['rho']*np.heaviside(self.A_grid_big,0.5))
        
    def run(self,mask_input_A,mask_input_R):
        tic = time.time()
        assert mask_input_A.shape[0] == self.Tsteps, 'input has incorrect length'
        assert mask_input_R.shape[0] == self.Tsteps, 'input has incorrect length'
        j = 0
        for i in range(self.Tsteps):
            self.t += self.dt
            
            self.update(mask_input_A[i],mask_input_R[i])
            
            if ((i % (self.Tsteps//10)) == 0) and self.p_r == True:
                print((100*i/self.Tsteps),'% done')
                print('Time elapsed: ',time.time()-tic)
            
            if (i % (np.rint(self.save_every/self.dt).astype(int))) == 0:
                self.A_saved[:,j] = self.getAGrid()
                self.R_saved[:,j] = self.getRGrid()
                j += 1
    def getAGrid(self):
        
        return accumulate_arr((np.zeros(self.coords.size,dtype = int),self.coords),self.A,(self.A_grid.shape[0],self.A_grid.shape[0]))[0,:]
    
    def getRGrid(self):
        
        return accumulate_arr((np.zeros(self.coords.size,dtype = int),self.coords),self.R,(self.A_grid.shape[0],self.A_grid.shape[0]))[0,:]