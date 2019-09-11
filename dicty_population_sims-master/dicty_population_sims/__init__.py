#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 14:35:27 2018

@author: goldalex
"""

from __future__ import division
import scipy.sparse as sparse
import numpy as np
class Grid:
    def __init__(self,gridType = "Square", Shape = (1,1), edgeLength = (1,1)):
        #Build a grid of the specified type (square, rectangular, disordered) with the specified size (e.g. # of points) and edge length
        self.shape = Shape
        self.edgeL = edgeLength
        if gridType == "Square":
            self.grid=np.ones(self.shape)
    
            self.positions=np.ndarray.tolist(self.grid)
            for i in range(self.shape[0]):
                for j in range(self.shape[1]):
                    self.positions[i][j]=(i,j)
    
    def getGradientOp(self,dir,BCs):
        #Return an operator that calculates the gradient along a specified direction, with BCs determining what to do at the boundary
        return 0
    
    def getLaplacianOp(self,BCs = "Dirichlet"):
        #Return an operator that calculates the Laplacian given BCs, eg a stencil to convolve with a grid or just a matrix to multiply by
        #start by creating the basic tridiagonal matrix
        assert self.shape[0]==self.shape[1], "non-square grids not supported yet"
        assert self.edgeL[0]==self.edgeL[1], "non-square grids not supported yet"
        L=self.shape[0]
        h=self.edgeL[0]
        diag=np.ones(L**2)
        odiag=([0]+[1]*(L-1))*(L)
        out = np.diag(-4*diag,0)+np.diag(odiag[1:],-1)+np.diag(odiag[1:],1)+np.diag(diag[L:],L)+np.diag(diag[L:],-L)
        out = out/(h**2)
        return out

class Agent:
    def __init__(self,pos,state,agentParams):
        #initialize the agent with some position vector "pos" and some state vector "state"
        self.pos=pos
        self.state=state
        self.params=agentParams
    def update(self,signals,dt):
        #update the position and internal state of the agent, given some update rule TBD
        #Set up placeholder variables for the state increment and the position increment
        stateInc=np.zeros(np.shape(self.state))
        posInc=(0,0)
        
        #set the next state
        stateInc[0]=self.state[1]*dt
        stateInc[1]=-self.params['omega']**2*self.state[0]*dt
        
        #loop over the signals and increment them however needed
        for i in range(np.size(signals)):
            pass
        
        #increment the position and set the next state
        self.pos+=posInc
        self.state+=stateInc

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

