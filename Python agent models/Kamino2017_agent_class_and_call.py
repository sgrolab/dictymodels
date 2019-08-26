    # -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 12:35:31 2019

Kamino 2017 model A"

@author: Chuqiao Huyan

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
        
    
    def update(self, signals, dt, overriding_sig = 'none'):
        tau=self.Param['tau']   
        n=self.Param['n'] 
        K=self.Param['K'] 
        kt=self.Param['kt'] 
        gamma=self.Param['gamma'] 
        delta=self.Param['delta'] 
        rho=self.Param['rho'] 
        
        dxdt=((self.z_now+delta)-self.x_now)/tau
        dydt=(self.z_now+delta)**n/((self.z_now+delta)**n+(K*self.x_now)**n)-self.y_now
        
        x_next=self.x_now+dxdt*dt
        y_next=self.y_now+dydt*dt
        if isinstance(overriding_sig, str):
            dzdt=rho*kt*self.y_now-gamma*(self.z_now-signals)
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
#%% Adaptive spike for individual cells
tau=1.5
n=2
K=4
kt=2
delta=0.01
gamma=3
rho= 0.01 # population density, doesn't matter for single cells
AgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
x0=0.01
y0=0.05
z0=0.005
Kamino_agent=Kamino2017_agent([x0,y0,z0],AgentParam)

x_trace=[x0]
y_trace=[y0]

dt=0.001
t_tot=25
t=list(np.arange(0,t_tot,dt))

stim_time_step=3000 # at this time step input is applied
constant_signal=0.3; # randomly set 2019/7/17
signal_trace=np.zeros(len(t))
signal_trace[stim_time_step:] = constant_signal


for i in range(len(t)-1):
    x_now=x_trace[i]
    y_now=y_trace[i]
    x_next,y_next,z_next= Kamino_agent.update(1, dt, signal_trace[i])
    x_trace.append(x_next)
    y_trace.append(y_next)
        
   
# Convert into np array
x_trace = np.array(x_trace) # vectorize p_trace
y_trace = np.array(y_trace)
t_plot = np.array(t)
#%% Figure single cell adaptive spike 
label_font_size = 25
trace_width = 6.0

fig5 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot,signal_trace, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
# ax1.set_xlabel('Time',fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot,y_trace, color='g',linewidth=trace_width)
ax2.set_ylabel(r'$cAMP_{cyto}$, y',fontsize=label_font_size)
ax2.yaxis.label.set_color('g')

ax3 = fig5.add_subplot(grid[1:, 0], sharex=ax2, frameon=False)
line2=ax3.plot(t_plot,x_trace, color='b', linewidth=trace_width) # right axis
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")
ax3.set_ylabel('Adaptation variable, x' ,fontsize=label_font_size)
ax3.yaxis.label.set_color('b')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax3.set_xlabel('Time',fontsize=label_font_size)

plt.show()
#%% Fig 5D
        
        
#%% Fig 5E
tau=1.5
n=2
K=4
kt=2
delta=0.01

gamma_space=np.logspace(0, 2.0, num=6)
rho_space=np.logspace(0, 2.0, num=6)

# Initialize oscillation phase matrix, based on z trace
PkWdthMean = np.zeros((len(gamma_space), len(rho_space))) # PkWdthMean- mean oscillation time
PkPrmMean = np.zeros((len(gamma_space), len(rho_space))) # PkPrmMean - mean oscillation peak prominence
OscOrNot = np.zeros((len(gamma_space), len(rho_space))) # OscOrNot - 1 or 0, oscillatory/ nonoscillatory


dt=0.0005 
t_tot=150
t=np.arange(0,t_tot,dt)
signal_trace=np.zeros(len(t))
#%%
for j in range(len(gamma_space)):
    gamma=gamma_space[j]
    for k in range(len(rho_space)):
        rho=rho_space[k]

        Param={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
        
        x0=0.01
        y0=0.08
        z0=0.01
        x_trace=[x0]
        y_trace=[y0]
        z_trace=[z0]
        
        test_agent=Kamino2017_agent([x0,y0,z0],Param)

        for i in range(len(t)-1):
            x_now=x_trace[i]
            y_now=y_trace[i]
            signal_now=signal_trace[i]
            
            x_next,y_next,z_next=test_agent.update(signal_now,dt)
            x_trace.append(x_next)
            y_trace.append(y_next)
            z_trace.append(z_next)
            
#        # check simulation traces
#        t_plot = np.array(t)
#        plt.figure()
#        plt.plot(t_plot,x_trace,t_plot,y_trace,t_plot,z_trace)
#        plt.xlabel('Time')
#        plt.ylabel('x,y,z')
#        plt.title('Fig5D with gamma= '+str(gamma)+' rho= '+str(rho))
#        plt.gca().legend(('x','y','z'))
#        plt.show()
        
        z_trace=np.array(z_trace) # convert list to array
        z_trace_later=z_trace[math.floor(len(t)*0.5):] # the later part of trace
        PkPos, PkProperties = find_peaks(z_trace_later, prominence=(0.1,1000))
        # Check find_peaks
        # plt.plot(z_trace_later)
        # plt.plot(peaks, z_trace_later[peaks], "x")
        if PkPos.size: # if there is oscillation
            PkWdthMean[k,j]=dt*np.mean(np.diff(PkPos))
            PkPrmMean[k,j]=np.mean(PkProperties["prominences"])
            OscOrNot[k,j]=1


#%%
fig, ax = plt.subplots()
p = ax.pcolor(np.log10(gamma_space), np.log10(rho_space), OscOrNot, cmap='RdBu', vmin=abs(OscOrNot).min(), vmax=abs(OscOrNot).max())
ax.set_xlabel(r'Dilution rate $\gamma$, log scale')
ax.set_ylabel(r'Cell density $\rho$, log scale')
fig.suptitle('Oscillation or not')
cb = fig.colorbar(p)
#%%
PkWdthMean[PkWdthMean == 0] = 'nan'
MeanFR=1/PkWdthMean # mean oscillation frequency/ firing rate
fig, ax = plt.subplots()
p = ax.pcolor(np.log10(gamma_space), np.log10(rho_space), MeanFR, cmap='jet', vmin=np.nanmin(MeanFR), vmax=np.nanmax(MeanFR))
ax.set_xlabel(r'Dilution rate $\gamma$, log scale')
ax.set_ylabel(r'Cell density $\rho$, log scale')
fig.suptitle('Mean oscillation frequency/ firing rate')
cb = fig.colorbar(p)
#%%
PkPrmMean[PkPrmMean == 0] = 'nan'
fig, ax = plt.subplots()
p = ax.pcolor(np.log10(gamma_space), np.log10(rho_space), PkPrmMean, cmap='jet', vmin=np.nanmin(PkPrmMean), vmax=np.nanmax(PkPrmMean))
ax.set_xlabel(r'Dilution rate $\gamma$, log scale')
ax.set_ylabel(r'Cell density $\rho$, log scale')
fig.suptitle('Mean oscillation peak prominence')
cb = fig.colorbar(p)

#%% Fig 5F
tau=1.5
n=2
K=4
kt=2
delta=0
gamma=3
rho=1

dt=0.0005 
t_tot=140
t=np.arange(0,t_tot,dt)
signal_trace=np.zeros(len(t))
temp_array=signal_trace[int(round(50/dt)):]
z0=0.02
signal_trace[int(round(50/dt)):]=z0

Param={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
        
x0=0.01
y0=0.08
z0=0.01
x_trace=[x0]
y_trace=[y0]
z_trace=[z0]

test_agent=Kamino2017_agent([x0,y0,z0],Param)

for i in range(len(t)-1):
    x_now=x_trace[i]
    y_now=y_trace[i]
    signal_now=signal_trace[i]
    
    x_next,y_next,z_next=test_agent.update(signal_now,dt)
    x_trace.append(x_next)
    y_trace.append(y_next)
    z_trace.append(z_next)

# check simulation traces
t_plot = np.array(t)
plt.plot(t_plot,x_trace,t_plot,y_trace,t_plot,z_trace)
plt.axvline(x=50,linestyle='--')
plt.xlabel('Time')
plt.ylabel('x,y,z')
plt.title('Fig5F with z0= '+str(z0)+r', $\rho$= '+str(rho))
plt.gca().legend(('x','y','z'))
plt.show()


#%% Fig 5G
tau=1.5
n=2
K=4
kt=2
delta=0.01
gamma=3

z0_space=np.logspace(-2.0, 0, num=20)
rho_space=np.logspace(0, 2.0, num=20)

# Initialize oscillation phase matrix, based on z trace
PkWdthMean = np.zeros((len(z0_space), len(rho_space))) # PkWdthMean- mean oscillation time
PkPrmMean = np.zeros((len(z0_space), len(rho_space))) # PkPrmMean - mean oscillation peak prominence
OscOrNot = np.zeros((len(z0_space), len(rho_space))) # OscOrNot - 1 or 0, oscillatory/ nonoscillatory

dt=0.0005 
t_tot=150
t=np.arange(0,t_tot,dt)
signal_trace=np.zeros(len(t))
#%%
for j in range(len(z0_space)):
    z0=z0_space[j]
    temp_array=signal_trace[int(round(50/dt)):]
    signal_trace[int(round(50/dt)):]=z0
    
    for k in range(len(rho_space)):
        rho=rho_space[k]

        Param={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
        
        x0=0.01
        y0=0.08
        z0=0.01
        x_trace=[x0]
        y_trace=[y0]
        z_trace=[z0]
        
        test_agent=Kamino2017_population([x0,y0,z0],Param)

        for i in range(len(t)-1):
            x_now=x_trace[i]
            y_now=y_trace[i]
            signal_now=signal_trace[i]
            
            x_next,y_next,z_next=test_agent.update(signal_now,dt)
            x_trace.append(x_next)
            y_trace.append(y_next)
            z_trace.append(z_next)
        # check simulation traces
        #t_plot = np.array(t)
        #plt.plot(t_plot,x_trace,t_plot,y_trace,t_plot,z_trace)
        #plt.xlabel('Time')
        #plt.ylabel('x,y,z')
        #plt.title('Fig5D with gamma= '+str(gamma)+' rho= '+str(rho))
        #plt.gca().legend(('x','y','z'))
        #plt.show()
        z_trace=np.array(z_trace) # convert list to array
        z_trace_later=z_trace[math.floor(len(t)*0.5):] # the later part of trace
        PkPos, PkProperties = find_peaks(z_trace_later, prominence=(0.1,1000))
        # Check find_peaks
        # plt.plot(z_trace_later)
        # plt.plot(peaks, z_trace_later[peaks], "x")
        if PkPos.size: # if there is oscillation
            PkWdthMean[k,j]=dt*np.mean(np.diff(PkPos))
            PkPrmMean[k,j]=np.mean(PkProperties["prominences"])
            OscOrNot[k,j]=1


#%%
fig, ax = plt.subplots()
p = ax.pcolor(np.log10(z0_space), np.log10(rho_space), OscOrNot, cmap='RdBu', vmin=abs(OscOrNot).min(), vmax=abs(OscOrNot).max())
ax.set_xlabel('Extracellular cAMP z0, log scale')
ax.set_ylabel(r'Cell density $\rho$, log scale')
fig.suptitle('Oscillation or not')
cb = fig.colorbar(p)
#%%
PkWdthMean[PkWdthMean == 0] = 'nan'
MeanFR=1/PkWdthMean # mean oscillation frequency/ firing rate
fig, ax = plt.subplots()
p = ax.pcolor(np.log10(z0_space), np.log10(rho_space), MeanFR, cmap='jet', vmin=np.nanmin(MeanFR), vmax=np.nanmax(MeanFR))
ax.set_xlabel('Extracellular cAMP z0, log scale')
ax.set_ylabel(r'Cell density $\rho$, log scale')
fig.suptitle('Mean oscillation frequency/ firing rate')
cb = fig.colorbar(p)
#%%
PkPrmMean[PkPrmMean == 0] = 'nan'
fig, ax = plt.subplots()
p = ax.pcolor(np.log10(z0_space), np.log10(rho_space), PkPrmMean, cmap='jet', vmin=np.nanmin(PkPrmMean), vmax=np.nanmax(PkPrmMean))
ax.set_xlabel('Extracellular cAMP z0, log scale')
ax.set_ylabel(r'Cell density $\rho$, log scale')
fig.suptitle('Mean oscillation peak prominence')
cb = fig.colorbar(p)
