#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 14:53:43 2019

@author: goldalex
"""

from __future__ import division
import numpy as np
import scipy as sp
import time
import matplotlib.pyplot as plt
from activeGelFramework import calc_square_laplacian_noflux_matrix
from scipy.ndimage import zoom
import matplotlib.animation as animation

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
        'box_size' : 50,
        'agent_dim' : .5,
        'num_agents' : 100*100,
        }

def scale(shape, B, k):     # scale B to shape scaled by k
    Y = shape[0]
    X = shape[1]
    A = np.zeros(shape)
    for y in range(0, k):
        for x in range(0, k):
            A[y:Y:k, x:X:k] = B
    
#    #cut off last row and column, we dont use them
#    A = A[0:-1,0:-1]
#    
#    avg_mat_diag = np.array(([1,0]*(size[1]//2)+[1]+[0]*size[1])*(size[0]//2)+[1])
#    avg_mat_up_one = np.array(([0,1/2]*(size[1]//2)+[0]+[0]*size[1])*(size[0]//2)+[0,1/2]*(size[1]//2))
#    avg_mat_down_one = np.array(([1/2,0]*(size[1]//2)+[0]+[0]*size[1])*(size[0]//2)+[1/2,0]*(size[1]//2))
#    
#    for i in range(shape[0]//2):
#        
#    
    return A

def run_FHN_2Grid_simulation(T,save_every,g_params=grid_params_default,c_params=cell_params_default,init_signal=False,init_cell_state=False, noise = False):
    
    dt=g_params['dx']**2/(8*g_params['D_sig'])
    Tsteps=np.int(T/dt)
    signal_size=np.shape(np.eye(int(g_params['box_size']/g_params['dx'])))
    agent_ratio = int(g_params['agent_dim']/g_params['dx'])
    agent_size=(int(g_params['box_size']/g_params['agent_dim']),int(g_params['box_size']/g_params['agent_dim']))
    
    if init_signal is False:
        signal = np.zeros(signal_size)
    else:
        if init_signal.shape == signal_size:
            signal = np.copy(init_signal)
        else:
            assert False, "Initial signal is the wrong shape"
    
    if init_cell_state is False:
        #initialize fields
        cell_state = np.ones((agent_size[0],agent_size[1],2))
        
        #define random block initial conditions
        ic = np.random.normal(loc = 0.0, scale = 2, size =(2,2,2))
        
        #implement random initial conditions in 2x2 coarse grained block, for activator and repressor
        cell_state[0:agent_size[0]//2,0:agent_size[1]//2] = ic[0,0]
        cell_state[agent_size[0]//2:,0:agent_size[1]//2] = ic[1,0]
        cell_state[0:agent_size[0]//2,agent_size[1]//2:] = ic[0,1]
        cell_state[agent_size[0]//2:,agent_size[1]//2:] = ic[1,1]
        A = cell_state[:,:,0]
        R = cell_state[:,:,1]
        
    else: 
        if init_cell_state.shape == (agent_size[0],agent_size[1],2):
            A = np.copy(init_cell_state[:,:,0])
            R = np.copy(init_cell_state[:,:,1])
        else:
            assert False, "Initial cell state vector has the wrong shape"
            
    #get laplacian matrix
    lap_mat = calc_square_laplacian_noflux_matrix(signal_size)
    
    #get interpolation matrix
    interp_mat = scale_down_mat_sp(signal_size,agent_ratio)
    
    #define probability to have an agent grid point to be vacant
    p_vacant = g_params['p_vacant']
    
    #define agent vacancy mask
    vacancy = np.random.choice(2,size=agent_size,p=[p_vacant,1-p_vacant])
    A = np.multiply(vacancy,A)
    R = np.multiply(vacancy,R)
    
    #save values
    A_saved = np.zeros((A.shape[0],A.shape[1],np.int(T/save_every)))
    R_saved = np.zeros((R.shape[0],R.shape[1],np.int(T/save_every)))
    signal_saved = np.zeros((signal.shape[0],signal.shape[1],np.int(T/save_every)))
    
    t=0
    j=0
    #begin timesteping
    for i in range(Tsteps):
        t += dt
        #prepare to calculate laplacian of cext

        #convolve to calculate laplacian
        laplacian = 1/(g_params['dx']**2)*np.reshape(lap_mat.dot(np.reshape(signal,signal_size[0]*signal_size[1])),(signal_size[0],signal_size[1]))
        
        #update
        #set up dummy update variables
        deltaCext = dt*(g_params['D_sig']*laplacian - c_params['J']*signal + c_params['rho']*c_params['a0'])\
        + dt*(np.multiply(np.repeat(np.repeat(vacancy,agent_ratio,0),agent_ratio,1),c_params['D']*c_params['rho']*np.heaviside(np.repeat(np.repeat(A,agent_ratio,0),agent_ratio,1),0.5)))
        
        deltaA = ((A-(np.power(A,3))/3 - R ))*dt\
        + ((c_params['a']*np.log1p(np.reshape(interp_mat.dot(np.reshape(signal,signal_size[0]*signal_size[1])),agent_size)/c_params['Kd'])))*dt 
        
        if noise is True:
            deltaA += + c_params['sigma']*np.random.normal(loc = 0.0, scale = np.sqrt(dt), size = agent_size)
        
        deltaR = (((A-c_params['gamma']*R)+c_params['c0'])*c_params['epsilon'])*dt
        
        #update for real
        signal += deltaCext
        A += np.multiply(deltaA,vacancy)
        R += np.multiply(deltaR,vacancy)
        
        if t >= save_every:
            A_saved[:,:,j] = A
            R_saved[:,:,j] = R
            signal_saved[:,:,j] = signal

            j += 1
            t = 0
        if i % (Tsteps//10) == 0:
            print(100*i/Tsteps,"% done")
        
    return (A_saved, R_saved, signal_saved)

def run_FHN_2Grid_simulation_diagnostics(T,save_every,g_params=grid_params_default,c_params=cell_params_default,init_signal=False,init_cell_state=False, noise = False):
    
    dt=g_params['dx']**2/(8*g_params['D_sig'])
    Tsteps=np.int(T/dt)
    signal_size=np.shape(np.eye(int(g_params['box_size']/g_params['dx'])))
    agent_ratio = int(g_params['agent_dim']/g_params['dx'])
    agent_size=(int(g_params['box_size']/g_params['agent_dim']),int(g_params['box_size']/g_params['agent_dim']))
    
    #set up diagnostics
    diagnostic_labels = ("time spend calculating laplacian","time spent updating signal","time spend calculating autonomous activator dynamics",\
                         "time spent calculating signal inputs to activator","time spend calculating activator noise", "time spent on repressor dynamics","time spend saving data")
    diagnostics = np.zeros((7))
    
    if init_signal is False:
        signal = np.zeros(signal_size)
    else:
        if init_signal.shape == signal_size:
            signal = np.copy(init_signal)
        else:
            assert False, "Initial signal is the wrong shape"
    
    if init_cell_state is False:
        #initialize fields
        cell_state = np.ones((agent_size[0],agent_size[1],2))
        
        #define random block initial conditions
        ic = np.random.normal(loc = 0.0, scale = 2, size =(2,2,2))
        
        #implement random initial conditions in 2x2 coarse grained block, for activator and repressor
        cell_state[0:agent_size[0]//2,0:agent_size[1]//2] = ic[0,0]
        cell_state[agent_size[0]//2:,0:agent_size[1]//2] = ic[1,0]
        cell_state[0:agent_size[0]//2,agent_size[1]//2:] = ic[0,1]
        cell_state[agent_size[0]//2:,agent_size[1]//2:] = ic[1,1]
        A = cell_state[:,:,0]
        R = cell_state[:,:,1]
        
    else: 
        if init_cell_state.shape == (agent_size[0],agent_size[1],2):
            A = np.copy(init_cell_state[:,:,0])
            R = np.copy(init_cell_state[:,:,1])
        else:
            assert False, "Initial cell state vector has the wrong shape"
            
    #get laplacian matrix
    lap_mat = calc_square_laplacian_noflux_matrix(signal_size)
    
    #define probability to have an agent grid point to be vacant
    p_vacant = g_params['p_vacant']
    
    #define agent vacancy mask
    vacancy = np.random.choice(2,size=agent_size,p=[p_vacant,1-p_vacant])
    A = np.multiply(vacancy,A)
    R = np.multiply(vacancy,R)
    
    #save values
    A_saved = np.expand_dims(np.array(A),2)
    R_saved = np.expand_dims(np.array(R),2)
    signal_saved = np.expand_dims(np.copy(signal),2)
    
    t=0
    j=0
    #begin timesteping
    for i in range(Tsteps):
        t += dt
        #prepare to calculate laplacian of cext
        lap_tic = time.time()
        #convolve to calculate laplacian
        laplacian = 1/(g_params['dx']**2)*np.reshape(lap_mat.dot(np.reshape(signal,signal_size[0]*signal_size[1])),(signal_size[0],signal_size[1]))
        diagnostics[0] += time.time() - lap_tic
        #update
        #set up dummy update variables
        signal_tic = time.time()
        deltaCext = dt*(g_params['D_sig']*laplacian - c_params['J']*signal + c_params['rho']*c_params['a0']) + dt*(np.multiply(zoom(vacancy,agent_ratio,order=1),c_params['D']*c_params['rho']*np.heaviside(zoom(A,agent_ratio,order=1),0.5)))
        diagnostics[1] += time.time() - signal_tic
        
        A_aut_tic = time.time()
        deltaA = ((A-(A*A*A)/3 - R ))*dt 
        diagnostics[2] += time.time() - A_aut_tic
        
        A_sig_tic = time.time()
        deltaA += ((c_params['a']*np.log1p(signal[0:signal_size[0]:agent_ratio,0:signal_size[1]:agent_ratio]/c_params['Kd'])))*dt 
        diagnostics[3] += time.time() - A_sig_tic
        
        A_noise_tic = time.time()
        if noise is True:
            deltaA += + c_params['sigma']*np.random.normal(loc = 0.0, scale = np.sqrt(dt), size = agent_size)
        diagnostics[4] += time.time() - A_noise_tic
        
        rep_tic = time.time()
        deltaR = (((A-c_params['gamma']*R)+c_params['c0'])*c_params['epsilon'])*dt
        diagnostics[5] += time.time() - rep_tic
        
        #update for real
        signal += deltaCext
        A += np.multiply(deltaA,vacancy)
        R += np.multiply(deltaR,vacancy)
        
        save_tic = time.time()
        if True:
            j += 1
            t = 0
            A_saved = np.append(A_saved,np.expand_dims(A,2),2)
            R_saved = np.append(R_saved,np.expand_dims(R,2),2)
            signal_saved = np.append(signal_saved,np.expand_dims(signal,2),2)
        if i % (Tsteps/10) == 0:
            print(100*i/Tsteps,"% done")
            
        diagnostics[6] += time.time()-save_tic
        
    return (A_saved, R_saved, signal_saved,diagnostic_labels,diagnostics)


def run_FHN_agents_masked(T,save_every,g_params=grid_params_default_masked,c_params=cell_params_default,init_signal=False,init_cell_state=False,init_coords = False, noise = False):
    
    dt=g_params['dx']**2/(8*g_params['D_sig'])
    Tsteps=np.int(T/dt)
    signal_size=np.shape(np.eye(int(g_params['box_size']/g_params['dx'])))
    agent_ratio = int(g_params['agent_dim']/g_params['dx'])
    
    if init_signal is False:
        signal = np.zeros(signal_size)
    else:
        if init_signal.shape == signal_size:
            signal = np.copy(init_signal)
        else:
            assert False, "Initial signal is the wrong shape"
    
    if init_cell_state is False:
        #initialize fields
        cell_state = np.random.normal(loc = 0.0, scale = 2, size = (2,g_params['num_agents']))
        
        A = cell_state[0,:]
        R = cell_state[1,:]
        
    else: 
        if init_cell_state.shape == (2,g_params['num_agents']):
            A = np.copy(init_cell_state[0,:])
            R = np.copy(init_cell_state[1,:])
        else:
            assert False, "Initial cell state vector has the wrong shape"
            
    #define agent locations
    agent_size=(int(g_params['box_size']/g_params['agent_dim']),int(g_params['box_size']/g_params['agent_dim']))
    
    if init_coords is False:
        coords = np.zeros((2,g_params['num_agents']),dtype = 'int')
        for agent in range(g_params['num_agents']):
            coords[:,agent] = [agent//agent_size[0],agent%agent_size[0]]
    else:
        assert init_coords.shape[1] == g_params['num_agents'], "numbers of coordinates and agents do not match"
        coords = init_coords.astype(int)
    
    A_grid = np.zeros((signal_size[0]//agent_ratio,signal_size[1]//agent_ratio))
    #get laplacian matrix
    lap_mat = calc_square_laplacian_noflux_matrix(signal_size)
    
    #get signal-to-agent scaledown matrix
    interp_mat = scale_down_mat_sp(signal_size,agent_ratio)
    
    #save values
    A_saved = np.zeros((A.shape[0],np.int(T/save_every)))
    R_saved = np.zeros((R.shape[0],np.int(T/save_every)))
    signal_saved = np.zeros((signal.shape[0],signal.shape[1],np.int(T/save_every)))
    coords_saved = np.zeros((coords.shape[0],coords.shape[1],np.int(T/save_every)),dtype = int)
    
    t=0
    j=0
    #begin timesteping
    for i in range(Tsteps):
        t += dt
        #prepare to calculate laplacian of cext

        #laplacian of signal
        laplacian = 1/(g_params['dx']**2)*np.reshape(lap_mat.dot(np.reshape(signal,signal_size[0]*signal_size[1])),(signal_size[0],signal_size[1]))
        
        #get activator grid
        try:
            A_grid = accumulate_arr(coords,A,A_grid.shape)
        except(ValueError):
            print(coords.shape)
            print(A.shape)
            assert False
        
        #update
        #set up dummy update variables
        deltaCext = dt*(g_params['D_sig']*laplacian - c_params['J']*signal + c_params['rho']*c_params['a0'])\
        + dt*(c_params['D']*c_params['rho']*np.heaviside(np.repeat(np.repeat(A_grid,agent_ratio,0),agent_ratio,1),0.5))
        
        deltaA = ((A-(np.power(A,3))/3 - R ))*dt\
        + ((c_params['a']*np.log1p(np.reshape(interp_mat.dot(np.reshape(signal,signal_size[0]*signal_size[1])),agent_size)[coords[0,:],coords[1,:]]/c_params['Kd'])))*dt 
        
        if noise is True:
            deltaA += + c_params['sigma']*np.random.normal(loc = 0.0, scale = np.sqrt(dt), size = A.shape)
        
        deltaR = (((A-c_params['gamma']*R)+c_params['c0'])*c_params['epsilon'])*dt
        
        #update for real
        signal += deltaCext
        A += deltaA
        R += deltaR
        
        if t >= save_every:
            A_saved[:,j] = A
            R_saved[:,j] = R
            signal_saved[:,:,j] = signal
            coords_saved[:,:,j] = coords
            
            j += 1
            t = 0
        if i % (Tsteps//10) == 0:
            print(100*i/Tsteps,"% done")
    return (A_saved, R_saved, signal_saved,coords_saved)

def accumulate_arr(coords, arr, shape):

    # Get linear indices to be used as IDs with bincount
    lidx = np.ravel_multi_index(coords, shape)
    # Or lidx = coords[0]*(coords[1].max()+1) + coords[1]

    # Accumulate arr with IDs from lidx
    return np.bincount(lidx,arr,minlength=shape[0]*shape[1]).reshape(shape[0],shape[1])

def save_stack(stack,path,vmax=2.1,vmin=-2.9):
    num_imgs = stack.shape[2]
    for i in range(num_imgs):
        fig, ax = plt.subplots()
        img1 = ax.imshow(stack[:,:,i], vmin=vmin, vmax=vmax)
        fig.colorbar(img1, ax=ax)
        p = "picsOut/"+ path + "/image"+f"{i:05d}"+".png"
        plt.savefig(p)
        plt.close(fig)
        
def generate_phase_map(params):
    #generate phase variable map for a given set of single cell parameters
    
    #set up simulation details
    dt = 10**(-5)
    T = 200
    Tsteps = T//dt
    #set phase = 0 at the very beginning of a spike
    A = -1.0
    R = 2./3
    
    A_trace = np.zeros(Tsteps)
    R_trace = np.zeros(Tsteps)
    
    #do the integration
    
    for i in range(Tsteps):
        A_trace[0] = A
        R_trace[0] = R
        
        deltaA = dt*(A-np.power(A,3)/3 - R)
        deltaR = dt*(A-params['gamma']*R+params['c0'])*params['epsilon']
        
        A += deltaA
        R += deltaR
    
    #extract a single period
    
def zoom_mat(size,z):
    #construct top matrix block row by row
    out = sp.sparse.bsr_matrix(np.array(([1]+[0]*(size[1]*size[0]-1))))
    for i in range(z-1):
        next_row = sp.sparse.bsr_matrix(np.array(([1.-(i+1.)/z]+[(i+1.)/z]+[0]*(size[1]*size[0]-2))))
        out=sp.sparse.vstack((out,next_row))
    block = out.todense()
    for i in range(size[0]*size[1]-2):
        nextblock = sp.sparse.bsr_matrix(np.roll(block,i+1,axis=1))
        out = sp.sparse.vstack((out,nextblock))
    lastrow = np.zeros((size[0]*size[1]))
    lastrow[-1] = 1
    out = sp.sparse.vstack((out,sp.sparse.bsr_matrix(lastrow)))
    return out
    
def get_phase_data(A_saved,R_saved, c_saved, params):
    R_shifted = R_saved-((params['a']*np.log1p(c_saved/params['Kd'])))
    theta = np.arctan2(R_shifted,A_saved)
    theta_grad = np.gradient(theta,axis = (0,1))

    theta_grad_x = theta_grad[0]
    theta_grad_y = theta_grad[1]
    #remove jumps from gradient due to theta circling around
#    for i in range(theta_grad_x.shape[0]):
#        for j in range(theta_grad_x.shape[1]):
#            if np.abs(theta_grad_x[i,j]) >= 1:
#                theta_grad_x[i,j] += theta[i-1,j]
    with np.nditer(theta_grad_x, op_flags=['readwrite']) as it:
        for elem in it:
            if np.abs(elem) >= 2:
                elem[...] += -np.sign(elem)*np.pi
    
    return(theta,theta_grad_x,theta_grad_y)
    
def get_phase_data_1D(A_saved,R_saved, c_saved, params):
    R_shifted = R_saved-((params['a']*np.log1p(c_saved/params['Kd'])))
    theta = np.arctan2(R_shifted,A_saved)
    theta_grad = np.gradient(theta,axis = (0))

    #remove jumps from gradient due to theta circling around
#    for i in range(theta_grad_x.shape[0]):
#        for j in range(theta_grad_x.shape[1]):
#            if np.abs(theta_grad_x[i,j]) >= 1:
#                theta_grad_x[i,j] += theta[i-1,j]
    with np.nditer(theta_grad, op_flags=['readwrite']) as it:
        for elem in it:
            if np.abs(elem) >= 2:
                elem[...] += -np.sign(elem)*np.pi
    
    return(theta,theta_grad)
    
def plot_phase_data(theta,theta_grad):
    fig, ax = plt.subplots(nrows=1,ncols=2,sharex=True, sharey=True)
    ax[0].set_aspect(1)
    ax[1].set_aspect(1)
    img1 = ax[0].quiver(theta_grad_y[:,:,idx],theta_grad_x[:,:,idx])
    img2 = ax[1].imshow(theta[:,:,idx], vmin=-np.pi, vmax=np.pi)
    fig.colorbar(img2, ax=ax)
    
def scale_down_mat(current_shape,scale):
    assert current_shape[0]//scale == current_shape[0]/scale, "Shape not divisible by scale"
    assert current_shape[1]//scale == current_shape[1]/scale, "Shape not divisible by scale"
    
    out_shape = (current_shape[0]//scale,current_shape[1]//scale)
    
    line = np.array(([1/scale**2]*scale+[0]*(current_shape[0]-scale))*scale+[0]*(current_shape[0]*(current_shape[1]-scale)))
    out_mat = np.array(line)
    for i in range(out_shape[0]-1):
        newline = np.roll(line,scale*(i+1))
        out_mat = np.vstack((out_mat,newline))
    block = np.array(out_mat)
    for i in range(out_shape[1]-1):
        newblock = np.roll(block,scale*current_shape[0]*(i+1))
        out_mat=np.vstack((out_mat,newblock))
    return out_mat

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
    return out_mat

def run_FHN_2Grid_1D_simulation(T,save_every,g_params=grid_params_default,c_params=cell_params_default,init_signal=False,init_cell_state=False, noise = False):
    
    dt=g_params['dx']**2/(8*g_params['D_sig'])
    Tsteps=np.int(T/dt)
    signal_size=np.shape(np.ones(int(g_params['box_size']/g_params['dx'])))
    agent_ratio = int(g_params['agent_dim']/g_params['dx'])
    agent_size=(int(g_params['box_size']/g_params['agent_dim']))
    
    if init_signal is False:
        signal = np.zeros(signal_size)
    else:
        if init_signal.shape == signal_size:
            signal = np.copy(init_signal)
        else:
            assert False, "Initial signal is the wrong shape"
    
    if init_cell_state is False:
        #initialize fields
        cell_state = np.ones((agent_size,2))
        
        #define random block initial conditions
        ic = np.random.normal(loc = 0.0, scale = 2, size =(2,2))
        
        #implement random initial conditions in 2x2 coarse grained block, for activator and repressor
        cell_state[0:agent_size//2,0] = ic[0,0]
        cell_state[agent_size//2:,0] = ic[1,0]
        cell_state[0:agent_size//2,1] = ic[0,1]
        cell_state[agent_size//2:,1] = ic[1,1]
        A = cell_state[:,0]
        R = cell_state[:,1]
        
    else: 
        if init_cell_state.shape == (agent_size[0],agent_size[1],2):
            A = np.copy(init_cell_state[:,:,0])
            R = np.copy(init_cell_state[:,:,1])
        else:
            assert False, "Initial cell state vector has the wrong shape"
            
    #get laplacian matrix
    lap_mat = sp.sparse.diags([np.ones(signal_size[0]-1),-2*np.ones(signal_size[0]),np.ones(signal_size[0]-1)],[-1,0,1],format = 'csr')
    lap_mat += sp.sparse.diags(np.array([1]+[0]*(signal_size[0]-2)+[1]),format = 'csr')
    
    #get interpolation matrix
    line = np.zeros(signal_size)
    line[:agent_ratio] = 1/agent_ratio
    interp_mat = np.array(line)
    for i in range(agent_size-1):
        newline = np.roll(line,agent_ratio*(i+1))
        interp_mat = np.vstack((interp_mat,newline))
    
    #define probability to have an agent grid point to be vacant
    p_vacant = g_params['p_vacant']
    
    #define agent vacancy mask
    vacancy = np.random.choice(2,size=agent_size,p=[p_vacant,1-p_vacant])
    A = np.multiply(vacancy,A)
    R = np.multiply(vacancy,R)
    
    #save values
    A_saved = np.expand_dims(np.array(A),2)
    R_saved = np.expand_dims(np.array(R),2)
    signal_saved = np.expand_dims(np.copy(signal),2)
    
    t=0
    j=0
    k=0
    l=0
    #begin timesteping
    for i in range(Tsteps):
        t += dt
        #prepare to calculate laplacian of cext

        #convolve to calculate laplacian
        laplacian = 1/(g_params['dx']**2)*(lap_mat.dot(signal))
        
        #update
        #set up dummy update variables
        deltaCext = dt*(g_params['D_sig']*laplacian - c_params['J']*signal + c_params['rho']*c_params['a0'])\
        + dt*(np.multiply(np.repeat(vacancy,agent_ratio,0),c_params['D']*c_params['rho']*np.heaviside((np.repeat(A,agent_ratio)),0.5)))
        
        deltaA = ((A-(np.power(A,3))/3 - R ))*dt\
        + ((c_params['a']*np.log1p((interp_mat.dot(signal))/c_params['Kd'])))*dt 
        
        if noise is True:
            deltaA += + c_params['sigma']*np.random.normal(loc = 0.0, scale = np.sqrt(dt), size = agent_size)
        
        deltaR = (((A-c_params['gamma']*R)+c_params['c0'])*c_params['epsilon'])*dt
        
        #update for real
        signal += deltaCext
        A += np.multiply(deltaA,vacancy)
        R += np.multiply(deltaR,vacancy)
        k += 1
        if t >= save_every:
            j += 1
            t = 0
            A_saved = np.append(A_saved,np.expand_dims(A,1),1)
            R_saved = np.append(R_saved,np.expand_dims(R,1),1)
            signal_saved = np.append(signal_saved,np.expand_dims(signal,1),1)
        if k >= (Tsteps/(20)):
            l+=1
            print(5*l,"% done")
            k=0
        
    return (A_saved, R_saved, signal_saved)

def make_video_phase(imgs,path,name): 
 
    plt.rcParams['animation.ffmpeg_path'] = u'/Users/goldalex/anaconda3/bin/ffmpeg' # Add the path of ffmpeg here!!    
    n_images=imgs.shape[2] 
    assert (n_images>1)    
    
    fig, ax = plt.subplots()
    
    #lineR, = ax.plot(xaxis_data[0],R_data[0],'c-',label="resources")
    img = ax.imshow(imgs[:,:,0], vmin = -np.pi, vmax = np.pi, cmap = plt.cm.twilight, animated = True)    


    def updatefig(img_num): 
    
        #lineR.set_data(xaxis_data[img_num],R_data[img_num],'r-')
        
        img.set_data(imgs[:,:,img_num])

        return img
    
    
    ani = animation.FuncAnimation(fig, updatefig, np.arange(1, n_images), interval=300)
    
    ani.save("picsOut/"+ path + "/" + name + ".mp4")
    
    plt.close(fig)