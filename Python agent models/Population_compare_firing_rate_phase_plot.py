# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 20:08:16 2019

@author: Chuqiao

Population firing rate phase diagram (vs dilution rate & density)

"""
import pandas as pd
import numpy as np
import random
import math
import matplotlib.pyplot as plt

#%% Sgro 2015
# simulation results from matlab
pop_rate_path = r'E:\bu\trainings and conferences\qbio2019\Sgro_pop_rate.xlsx'
pop_rate = pd.read_excel(pop_rate_path)
pop_rate_mat_Sgro=pop_rate.as_matrix()

logrho_Sgro = np.linspace(-5.5,-3,num=26)
j_Sgro = np.linspace(0,1,num=21)

#%% Gregor 2010
pop_rate_path = r'E:\bu\trainings and conferences\qbio2019\Gregor_pop_rate.xlsx'
pop_rate_Gregor = pd.read_excel(pop_rate_path)
pop_rate_mat_Gregor=pop_rate_Gregor.as_matrix()

logrho_Gregor = np.linspace(-3.5,1,num=26)
k_Gregor = np.linspace(1,100,num=21)
#%% Kamino 2017
import numpy as np
import random
import math
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import find_peaks

# remove all variables in the work space
# sys.modules[__name__].__dict__.clear()


        
  #%%  
from Kamino2017_agent_and_pop_FUN import Kamino2017_pop
from time import perf_counter 
   
tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma_space=np.logspace(0, 2.0, num=3)
loggamma_space_Kamino=np.log10(gamma_space)
rho_space=np.logspace(0, 2.0, num=4)
logrho_space_Kamino = np.log10(rho_space)

# Initialize oscillation phase matrix, based on z trace
PkWdthMean = np.zeros((len(gamma_space), len(rho_space))) # PkWdthMean- mean oscillation time
PkPrmMean = np.zeros((len(gamma_space), len(rho_space))) # PkPrmMean - mean oscillation peak prominence

OscOrNot = np.zeros((len(gamma_space), len(rho_space))) # OscOrNot - 1 or 0, oscillatory/ nonoscillatory
pop_rate_Kamino = np.zeros((len(gamma_space), len(rho_space))) # population firing rate


dt=0.0005
t_tot=150
t=np.arange(0,t_tot,dt)
z0_influx = 0

# start simulation
tic = perf_counter() 
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
        
        test_agent=Kamino2017_pop([x0,y0,z0],Param)

        for i in range(len(t)-1):
            x_now=x_trace[i]
            y_now=y_trace[i]
            
            x_next,y_next,z_next=test_agent.update(z0_influx,dt)
            x_trace.append(x_next)
            y_trace.append(y_next)
            z_trace.append(z_next)
            
            #  check simulation traces
            fig = plt.figure()
            t_plot_Kamino = np.array(t)
            plt.plot(t_plot_Kamino,y_trace,t_plot_Kamino,z_trace)
            plt.xlabel('Time')
            plt.ylabel('x,y,z')
            plt.title('Fig5D with gamma= '+str(gamma)+' rho= '+str(rho))
            plt.gca().legend(('y','z'))
            plt.show()
        
        y_trace=np.array(y_trace) # convert list to array
        y_trace_later=y_trace[math.floor(len(t)*0.5):] # the later part of trace
        PkPos, PkProperties = find_peaks(y_trace_later, prominence=(0.02,100))
        # Check find_peaks
        fig = plt.figure()
        plt.plot(y_trace_later)
        plt.plot(PkPos, y_trace_later[PkPos], "x")
        if PkPos.size: # if there is oscillation
            # PkWdthMean[k,j]=dt*np.mean(np.diff(PkPos))
            # PkPrmMean[k,j]=np.mean(PkProperties["prominences"])
            OscOrNot[j,k]=1
            pop_rate_Kamino[j,k]=t_tot/len(PkProperties["prominences"])
    

    print('Finished gamma='+str(gamma))
    toc = perf_counter() 
    print('time passed'+str(toc-tic))
#%% 
title_font_size = 20 
label_font_size = 20
tick_font_size = 16

fig3 = plt.figure(figsize=(23, 4.5))
grid = plt.GridSpec(1,3, wspace=0.3, hspace=0.05)


ax2= fig3.add_subplot(grid[0,1])
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, OscOrNot.transpose(), cmap='jet') # cmap='jet'
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
heatmap = ax2.pcolor(gamma_space, logrho_space_Kamino, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
# ax2.set_xscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}(\gamma)$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Kamino & Sawai 2017', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
#%% Plot all
title_font_size = 20 
label_font_size = 20
tick_font_size = 16

fig3 = plt.figure(figsize=(23, 4.5))
grid = plt.GridSpec(1,3, wspace=0.3, hspace=0.05)

#ax1= fig3.add_subplot(grid[0,0])
#heatmap = ax1.pcolor(k_Gregor, logrho_Gregor, pop_rate_mat_Gregor, cmap='jet') # cmap='jet'
#heatmap.set_clim(0,0.16)
#cbar=fig3.colorbar(heatmap, ax=ax1);cbar.ax.tick_params(labelsize = tick_font_size) 
## ax1.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
## ax1.set_xlabel('Dilution rate J',fontsize=label_font_size)
#ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
#ax1.set_title('Gregor 2010', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax2= fig3.add_subplot(grid[0,1])
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, OscOrNot.transpose(), cmap='jet') # cmap='jet'
# heatmap = ax2.pcolor(loggamma_space_Kamino, logrho_space_Kamino, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
heatmap = ax2.pcolor(gamma_space, logrho_space_Kamino, pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
# ax2.set_xscale('log')
# heatmap.set_clim(0,0.16)
cbar=fig3.colorbar(heatmap, ax=ax2);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax2.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax2.set_xlabel(r'Dilution rate, $log_{10}(\gamma)$',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Kamino & Sawai 2017', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax3= fig3.add_subplot(grid[0,2])
heatmap = ax3.pcolor(100*j_Sgro, logrho_Sgro, pop_rate_mat_Sgro, cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.6)
cbar=fig3.colorbar(heatmap, ax=ax3);cbar.ax.tick_params(labelsize = tick_font_size) 
# ax3.set_ylabel(r'$log_{10}(\rho)$',fontsize=label_font_size)
# ax3.set_xlabel('Dilution rate',fontsize=label_font_size)
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Sgro 2015', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# fig3.text(0.5, 0.04,'Dilution rate', ha='center', va='center',fontsize=8)
fig3.text(0.06, 0.5, r'Cell density, $log_{10}(\rho)$', ha='center', va='center', rotation='vertical',fontsize=label_font_size)


plt.show()


#%% Kamino wired trend 
##### 7/29 7 pm need to run longer
from time import perf_counter 
   
tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma_space=np.logspace(0, 2.0, num=21)
loggamma_space_Kamino=np.log10(gamma_space)
rho_space=np.logspace(0, 2.0, num=26)
logrho_space_Kamino = np.log10(rho_space)

# Initialize oscillation phase matrix, based on z trace
PkWdthMean = np.zeros((len(gamma_space), len(rho_space))) # PkWdthMean- mean oscillation time
PkPrmMean = np.zeros((len(gamma_space), len(rho_space))) # PkPrmMean - mean oscillation peak prominence

OscOrNot = np.zeros((len(gamma_space), len(rho_space))) # OscOrNot - 1 or 0, oscillatory/ nonoscillatory
pop_rate_Kamino = np.zeros((len(gamma_space), len(rho_space))) # population firing rate


dt=0.0005 
t_tot=800
t=np.arange(0,t_tot,dt)
signal_trace=np.zeros(len(t)) # z0, background cAMP signal

j_test=[12]; k_test=[17,20];
# start si0mulation
tic = perf_counter() 
for j in range(len(gamma_space)): # j_test: #
    gamma=gamma_space[j]
    for k in range(len(rho_space)): # k_test: # 
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
            
#        #  check simulation traces
#        fig = plt.figure()
#        t_plot_Kamino = np.array(t)
#        plt.plot(t_plot_Kamino,y_trace,t_plot_Kamino,z_trace)
#        plt.xlabel('Time')
#        plt.ylabel('x,y,z')
#        plt.title('Fig5D with gamma= '+str(gamma)+' rho= '+str(rho))
#        plt.gca().legend(('y','z'))
#        plt.show()
        
        y_trace=np.array(y_trace) # convert list to array
        y_trace_later=y_trace[math.floor(len(t)*0.5):] # the later part of trace
        PkPos, PkProperties = find_peaks(y_trace_later, prominence=(0.02,100))
        
#         # Check find_peaks 
#        fig = plt.figure()
#        plt.plot(y_trace_later)
#        plt.plot(PkPos, y_trace_later[PkPos], "x")
#        plt.title('Population trace with gamma= '+str(gamma)+' rho= '+str(rho))
#        plt.show()
#         
        if PkPos.size: # if there is oscillation
            # PkWdthMean[k,j]=dt*np.mean(np.diff(PkPos))
            # PkPrmMean[k,j]=np.mean(PkProperties["prominences"])
            OscOrNot[j,k]=1
            pop_rate_Kamino[j,k]=len(PkProperties["prominences"])/t_tot
    

    print('Finished gamma='+str(gamma))
    toc = perf_counter() 
    print('time passed'+str(toc-tic))