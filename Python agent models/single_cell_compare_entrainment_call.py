# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:23:10 2019

@author: Chuqiao
"""
# entrainment quality
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import scipy.io

from scipy.signal import find_peaks, correlate 

# Normalization parameters
from NormParam import *
#%% Sgro 2015 with noise
from Sgro2015_agent_and_pop_FUN import Sgro2015_agent
# from TracesXcorr import TracesXcorr

#AdptPkPeriod=1.5
#period_space_Sgro =np.linspace(0.9/1.5*AdptPkPeriod,1.7/1.5*AdptPkPeriod,num=7) # period of cAMP stim 
#PkWdth_space_Sgro=  np.linspace(0.3/1.5*AdptPkPeriod, 1.4/1.5*AdptPkPeriod, num=8) # np.logspace(0.5, 0.5, num=1) 
#period_space_Sgro =np.linspace(0.5*AdptPkPeriod,1.15*AdptPkPeriod,num=7) # period of cAMP stim 
#PkWdth_space_Sgro=  np.linspace(0.2*AdptPkPeriod, 0.95*AdptPkPeriod, num=8) # np.logspace(0.5, 0.5, num=1) 

period_space_Sgro =np.linspace(0.8,1.8,num=8) # period of cAMP stim 
PkWdth_space_Sgro=  np.linspace(0.3, 1.4, num=10) # np.logspace(0.5, 0.5, num=1) 



# Initialize 
MeanR_Sgro = np.empty((len(PkWdth_space_Sgro), len(period_space_Sgro))) 
MeanR_Sgro[:] = np.nan # Mean cross-correlation of spike trace in cycle of stimulus, as nan matrix

MeanRnew_Sgro = np.empty((len(PkWdth_space_Sgro), len(period_space_Sgro))) 
MeanRnew_Sgro[:] = np.nan

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
A0=-1.5; R0=-0.5

dt=0.005
cAMP = 1 # extracellular cAMP
NumofCycle = 8 # cycle stimulation

j_test=[3]; k_test=[0,1,2,3,4]
for j in range(len(period_space_Sgro)): # j_test: #
    period = period_space_Sgro[j]
    for k in range(len(PkWdth_space_Sgro)): # k_test: #
        PkWdth = PkWdth_space_Sgro[k]
        if PkWdth>=period: 
            pass
        else:
            stim_1cycle=np.zeros(int(period*Nt_Sgro/dt))
            stim_1cycle[0:int(PkWdth*Nt_Sgro/dt)]=cAMP;
            signal_trace=np.concatenate((np.zeros(int(1*Nt_Sgro/dt)), np.tile(stim_1cycle,(NumofCycle)),np.zeros(int(1*Nt_Sgro/dt))),axis=0);
            
            t = list(np.arange(0,len(signal_trace)*dt,dt))
            t_plot_Sgro=np.array(t)/Nt_Sgro 
            
            A_trace_orig=[A0]; R_trace_orig=[R0]
            Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)
            for i in range(len(t)-1):
                A_now=A_trace_orig[i]
                R_now=R_trace_orig[i]
                signal_now=signal_trace[i]
        
                A_next,R_next,r_now=Sgro_agent.update(dt,signal_now)
                A_trace_orig.append(A_next)
                R_trace_orig.append(R_next)
       
            A_trace_offset=1.5
            A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
            A_trace_plot=(A_trace_orig+A_trace_offset)/Nh_Sgro;
            
            r=np.zeros(NumofCycle-1) # list that stores correlation coefficient to the first peak
            r_new = np.zeros(NumofCycle-1) # correlation coeff take in height difference
            InitPk=A_trace_plot[int((1+0*period)*Nt_Sgro/dt): int((1+1*period)*Nt_Sgro/dt)]

            for m in range(NumofCycle-1):
                FollowPk=A_trace_plot[int((1+(m+1)*period)*Nt_Sgro/dt) : (int((1+(m+1)*period)*Nt_Sgro/dt)+len(InitPk))]
                # SingleTraces= np.column_stack((SingleTraces,FollowPk)) 
                for m in range(NumofCycle-1):
                    FollowPk=A_trace_plot[int((1+(m+1)*period)*Nt_Sgro/dt) : (int((1+(m+1)*period)*Nt_Sgro/dt)+len(InitPk))]
                    R = np.corrcoef(InitPk,FollowPk)
                    r[m]=R[1,0]
                    r_new[m]=r[m]*(1-(sqrt(np.var(InitPk))-sqrt(np.var(FollowPk)))/sqrt(np.var(InitPk)))
            MeanR_Sgro[k,j] = np.mean(r)
            MeanRnew_Sgro[k,j] = np.mean(r_new)
            
#            # check signal_trace & response trace
#            fig3 = plt.figure(figsize=(4,4)); grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
#            ax1= fig3.add_subplot(grid[0, 0])
#            ax1.plot(t_plot_Sgro,signal_trace)
#            ax1.set_ylabel('extracellular cAMP'); 
#            ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#            ax2= fig3.add_subplot(grid[1:, 0])
#            ax2.plot(t_plot_Sgro,A_trace_plot)
#            fig3.text(0.5, 0.04, 'E-quality='+str(MeanR_Sgro[k,j])+', new quality=' 
#                      +str(MeanRnew_Sgro[k,j]),fontsize=10, ha='center')
#            plt.show()
            
#            # check InitPk vs.FollowPk
#            fig3 = plt.figure(figsize=(3,3)); grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
#            ax1= fig3.add_subplot(grid[0, 0])
#            ax1.plot(InitPk); ax1.plot(FollowPk)
#            ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#            plt.show()  
  
        
                    # MeanR_Sgro[k,j] = TracesXcorr(SingleTraces)

#% plotperiod- peak width- entrainment quality heatmap
label_font_size = 20
tick_font_size = 16
fig3 = plt.figure(figsize=(7, 6))
grid = plt.GridSpec(1, 1, wspace=0.1, hspace=0.1)

ax1= fig3.add_subplot(grid[0, 0])
# heatmap = ax1.pcolor(period_space_Sgro*10/np.amax(period_space_Sgro), PkWdth_space_Sgro* 1.4/1.7*10/np.amax(PkWdth_space_Sgro),MeanR_Sgro, cmap='jet') # cmap='jet'
heatmap = ax1.pcolor(period_space_Sgro, PkWdth_space_Sgro,MeanR_Sgro, cmap='jet')
heatmap.set_clim(0,1)
fig3.colorbar(heatmap, ax=ax1)

ax1.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Sgro2015 single cell NEW entrainment quality')
plt.show()


#%% Kamino 2017
from Kamino2017_agent_and_pop_FUN import Kamino2017_agent 

#AdptPkPeriod=8.5
#period_space_Kamino = np.linspace(0.5*AdptPkPeriod,1.15*AdptPkPeriod,num=7) # period of cAMP stim # np.array([5.4,7]) #
#PkWdth_space_Kamino=  np.linspace(0.2*AdptPkPeriod,0.95*AdptPkPeriod, num=8) # np.logspace(0.5, 0.5, num=1)  np.array([1.3,4]) # 
period_space_Kamino =np.linspace(0.8,1.8,num=8) # period of cAMP stim 
PkWdth_space_Kamino=  np.linspace(0.3, 1.4, num=10) # np.logspace(0.5, 0.5, num=1) 

# Initialize 
MeanR_Kamino = np.empty((len(PkWdth_space_Kamino), len(period_space_Kamino))) 
MeanR_Kamino[:] = np.nan # Mean cross-correlation of spike trace in cycle of stimulus, as nan matrix
MeanRnew_Kamino = np.empty((len(PkWdth_space_Kamino), len(period_space_Kamino))) 
MeanRnew_Kamino[:] = np.nan 

tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
KaminoAgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
x0=0.01; y0=0.05; z0=0.005

dt=0.001 ; # t_tot=20*Nt; t=list(np.arange(0,t_tot,dt))
cAMP = 1 # extracellular cAMP
NumofCycle = 8 # cycle stimulation

j_test=[1]; k_test=[0,1,2,3,4,5]
for j in range(len(period_space_Kamino)): # j_test: # 
    period = period_space_Kamino[j]
    for k in range(len(PkWdth_space_Kamino)): # k_test: # 
        PkWdth = PkWdth_space_Kamino[k]
        if PkWdth>=period: 
            pass
        else:
            stim_1cycle=np.zeros(int(period*Nt_Kamino/dt))
            stim_1cycle[0:int(PkWdth*Nt_Kamino/dt)]=cAMP;
            signal_trace=np.concatenate((np.zeros(int(1*Nt_Kamino/dt)), np.tile(stim_1cycle,(NumofCycle)),np.zeros(int(1*Nt_Kamino/dt))),axis=0);
            
            t = list(np.arange(0,len(signal_trace)*dt,dt))
            t_plot_Kamino=np.array(t)/Nt_Kamino
            
            x_trace=[x0]; y_trace=[y0]
            Kamino_agent=Kamino2017_agent([x0,y0,z0],KaminoAgentParam)
            for i in range(len(t)-1):
                x_now=x_trace[i]
                y_now=y_trace[i]
                x_next,y_next,z_next= Kamino_agent.update(dt, signal_trace[i])
                x_trace.append(x_next)
                y_trace.append(y_next)
                
       
            y_trace = np.array(y_trace)# vectorize y_trace
            y_trace = y_trace/Nh_Kamino

            r=np.zeros(NumofCycle-1) # list that stores correlation coefficient to the first peak
            r_new = np.zeros(NumofCycle-1) # correlation coeff take in height difference
            InitPk=y_trace[int(1*Nt_Kamino/dt): int((1+period)*Nt_Kamino/dt)]
            for m in range(NumofCycle-1):
                FollowPk=y_trace[int((1+(m+1)*period)*Nt_Kamino/dt) : (int((1+(m+1)*period)*Nt_Kamino/dt)+len(InitPk))]
                R = np.corrcoef(InitPk,FollowPk)
                r[m]=R[1,0]
                r_new[m]=r[m]*(1-(sqrt(np.var(InitPk))-sqrt(np.var(FollowPk)))/sqrt(np.var(InitPk)))
#                # check InitPk vs.FollowPk
#                fig3 = plt.figure(figsize=(3,3)); grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
#                ax1= fig3.add_subplot(grid[0, 0])
#                ax1.plot(InitPk); ax1.plot(FollowPk)
#                ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#                plt.show()
            MeanR_Kamino[k,j] = np.mean(r)
            MeanRnew_Kamino[k,j] = np.mean(r_new)

#             # check signal_trace & response trace
#            fig3 = plt.figure(figsize=(4,4)); grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
#            ax1= fig3.add_subplot(grid[0, 0])
#            ax1.plot(t_plot_Kamino,signal_trace)
#            ax1.set_ylabel('extracellular cAMP'); 
#            ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#            ax2= fig3.add_subplot(grid[1:, 0])
#            ax2.plot(t_plot_Kamino,y_trace)
#            fig3.text(0.5, 0.04, 'E-quality='+str(MeanR_Kamino[k,j])+', new quality=' 
#                      +str(MeanRnew_Kamino[k,j]),fontsize=10, ha='center')
#            plt.show()
            

#%plotperiod- peak width- entrainment quality heatmap
label_font_size = 20
tick_font_size = 16
# fig3 = plt.figure(figsize=(7, 6))
fig3 = plt.figure(figsize=(7, 6))
grid = plt.GridSpec(1, 1, wspace=0.1, hspace=0.1)

ax1= fig3.add_subplot(grid[0, 0])
heatmap = ax1.pcolor(period_space_Kamino, PkWdth_space_Kamino,MeanR_Kamino, cmap='jet')
heatmap.set_clim(0,1)
fig3.colorbar(heatmap, ax=ax1)
ax1.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Kamino2017 single cell NEW entrainment quality')
plt.show()


#%% Goldbeter 1987
from Goldbeter1987_agent_and_pop_FUN import Goldbeter1987_agent_3var

#AdptPkPeriod = 7.5
#period_space_Gold = np.linspace(0.5*AdptPkPeriod,1.15*AdptPkPeriod,num=7) # period of cAMP stim # np.array([5.4,7]) #
#PkWdth_space_Gold=  np.linspace(0.2*AdptPkPeriod,0.95*AdptPkPeriod, num=8) # np.logspace(0.5, 0.5, num=1)  np.array([1.3,4]) # 

period_space_Gold =np.linspace(0.8,1.8,num=8) # period of cAMP stim 
PkWdth_space_Gold = np.linspace(0.3, 1.4, num=10) # np.logspace(0.5, 0.5, num=1) 

# Initialize 
MeanR_Gold = np.empty((len(PkWdth_space_Gold), len(period_space_Gold))) 
MeanR_Gold[:] = np.nan # Mean cross-correlation of spike trace in cycle of stimulus, as nan matrix
MeanRnew_Gold = np.empty((len(PkWdth_space_Gold), len(period_space_Gold))) 
MeanRnew_Gold[:] = np.nan

k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10; L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01; theta=0.01
e= 1
q=4000
sig=0.6 # compared to 0.57
v=12; k= 4 # k prime in the paper
ki=1.7 
kt=0.9
kc=5.4 
h=5

Goldbeter3AgentParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}

p0=0.8; a0=3; b0=0.9; g0=0

dt=0.001 ; # t_tot=20*Nt; t=list(np.arange(0,t_tot,dt))
cAMP = 1 # extracellular cAMP
NumofCycle = 8 # cycle stimulation

#j_test=[0]; k_test=[0,3,5]
#period_space = [1];PkWdth_space = [0.3,0.6,0.9]
#count=0
fig3 = plt.figure(figsize=(8,8)); grid = plt.GridSpec(3, 3, wspace=0.3, hspace=0.3)
for  j in range(len(period_space_Gold)):# period in period_space: #j_test:# 
    period = period_space_Gold[j]
    for k in range(len(PkWdth_space_Gold)): # PkWdth in PkWdth_space: #k in k_test:#
        PkWdth = PkWdth_space_Gold[k]
        if PkWdth>=period: 
            pass
        else:
            stim_1cycle=np.zeros(int(period*Nt_Goldbeter/dt))
            stim_1cycle[0:int(PkWdth*Nt_Goldbeter/dt)]=cAMP;
            signal_trace=np.concatenate((np.zeros(int(1*Nt_Goldbeter/dt)), np.tile(stim_1cycle,(NumofCycle)),np.zeros(int(1*Nt_Goldbeter/dt))),axis=0);
            
            t = list(np.arange(0,len(signal_trace)*dt,dt))
            t_plot_Gold=np.array(t)/Nt_Goldbeter
            
            p_trace=[p0]; b_trace=[b0]; g_trace=[g0]
            Goldbeter3_agent=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
            for i in range(len(t)-1):
                p_now=p_trace[i]
                p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,signal_trace[i])
                p_trace.append(p_next)
                b_trace.append(b_next)
                g_trace.append(g_next)
                      
            b_trace = np.array(b_trace); b_trace = b_trace/Nh_Goldbeter

            r=np.zeros(NumofCycle-1) # list that stores correlation coefficient to the first peak
            r_new = np.zeros(NumofCycle-1) # correlation coeff take in height difference
            InitPk=b_trace[int(1*Nt_Goldbeter/dt): int((1+period)*Nt_Goldbeter/dt)]
            for m in range(NumofCycle-1):
                FollowPk=b_trace[int((1+(m+1)*period)*Nt_Goldbeter/dt) : (int((1+(m+1)*period)*Nt_Goldbeter/dt)+len(InitPk))]
                R = np.corrcoef(InitPk,FollowPk)
                r[m]=R[1,0]
                r_new[m]=r[m]*(1-(sqrt(np.var(InitPk))-sqrt(np.var(FollowPk)))/sqrt(np.var(InitPk)))
                    
            MeanR_Gold[k,j] = np.mean(r)
            MeanRnew_Gold[k,j] = np.mean(r_new)
            
#            # check signal_trace & response trace
#            fig3 = plt.figure(figsize=(5,5)); grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
#            ax1= fig3.add_subplot(grid[0, 0])
#            ax1.plot(t_plot_Gold,signal_trace)
#            ax1.set_ylabel('extracellular cAMP'); ax1.set_xlabel('Time, A.U.');
#            ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#            ax2= fig3.add_subplot(grid[1:, 0])
#            ax2.plot(t_plot_Gold,b_trace)
#            fig3.text(0.5, 0.04, 'E-quality='+str(MeanR_Gold[k,j])+', new quality=' 
#                      +str(MeanRnew_Gold[k,j]),fontsize=10, ha='center')
#            plt.show()      
#            # plot supplimentary figure         
#            ax1= fig3.add_subplot(grid[0, count])
#            ax1.plot(t_plot_Gold,signal_trace)
#            ax1.set_ylabel('extracellular cAMP'); ax1.set_xlabel('Time, A.U.');
#            ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#            ax2= fig3.add_subplot(grid[1:, count])
#            ax2.plot(t_plot_Gold,b_trace)
#            fig3.text(0.2+count*0.3, 0.02, 'Modified ntrainment quality=\n {:.4f}'.format(MeanRnew_Gold[k,j]),fontsize=12, ha='center')
#            plt.show() 
#            count=count+1
            #                # check InitPk vs.FollowPk
#                fig3 = plt.figure(figsize=(3,3)); grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
#                ax1= fig3.add_subplot(grid[0, 0])
#                ax1.plot(InitPk); ax1.plot(FollowPk)
#                ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#                plt.show()   
                
#% plotperiod- peak width- entrainment quality heatmap
label_font_size = 16
tick_font_size = 16
# fig3 = plt.figure(figsize=(7, 6))
fig3 = plt.figure(figsize=(7, 6))
grid = plt.GridSpec(1, 1, wspace=0.1, hspace=0.5)

ax1= fig3.add_subplot(grid[0, 0])
heatmap = ax1.pcolor(period_space_Gold, PkWdth_space_Gold,MeanR_Gold, cmap='jet')
heatmap.set_clim(0,1)
fig3.colorbar(heatmap, ax=ax1)
ax1.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Modified Entrainment period, A.U.',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Goldbeter 1987 single cell\n modified entrainment quality',fontsize = 16)
plt.show()

        
   
#%% Maeda and Loomis 2004
from MaedaLoomis2004_agent_and_pop_FUN import MaedaLoomis2004_agent


period_space_Maeda = np.linspace(0.8,1.8,num=8) # period of cAMP stim 
PkWdth_space_Maeda =  np.linspace(0.3, 1.4, num=10) 
# Initialize 
MeanR_Maeda = np.empty((len(PkWdth_space_Maeda), len(period_space_Maeda))) 
MeanR_Maeda[:] = np.nan # Mean cross-correlation of spike trace in cycle of stimulus, as nan matrix
MeanRnew_Maeda = np.empty((len(PkWdth_space_Maeda), len(period_space_Maeda))) 
MeanRnew_Maeda[:] = np.nan 

k1=2.0; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=1.0; k8=1.3; k9=0.3; k10=0.8
k11=0.7; k12=4.9; k13=23; k14=4.5
MaedaAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}
ACA0=0.1; PKA0=0.1; ERK20=0.1; RegA0=0.1; cAMPi0=0.01; 
cAMPe0=0.1; CAR10=0.1
state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]

dt=0.0005
cAMP = 1 # extracellular cAMP
NumofCycle = 8 # cycle stimulation

j_test=[0] ;k_test=[0,1,2,3,4,5]; 

for j in range(len(period_space_Maeda)):# j_test:# 
    period = period_space_Maeda[j]
    for k in  range(len(PkWdth_space_Maeda)): # k_test:#
        PkWdth = PkWdth_space_Maeda[k]
        if PkWdth>=period: 
            pass
        else:
            stim_1cycle=np.zeros(int(period*Nt_Maeda/dt))
            stim_1cycle[0:int(PkWdth*Nt_Maeda/dt)]=cAMP;
            signal_trace=np.concatenate((np.zeros(int(1*Nt_Maeda/dt)), np.tile(stim_1cycle,(NumofCycle)),np.zeros(int(1*Nt_Maeda/dt))),axis=0);
            
            t = list(np.arange(0,len(signal_trace)*dt,dt))
            t_plot_Maeda=np.array(t)/Nt_Maeda
            
            ACA_trace=[ACA0]; PKA_trace=[PKA0]; ERK2_trace=[ERK20]; RegA_trace=[RegA0]; 
            cAMPi_trace=[cAMPi0]; cAMPe_trace=[cAMPe0]; CAR1_trace=[CAR10]
            Maeda_agent=MaedaLoomis2004_agent([1,1],state0,MaedaAgentParam)
            for i in range(len(t)-1):
                ACA_now=ACA_trace[i]
                PKA_now=PKA_trace[i]
                ERK2_now=ERK2_trace[i]
                RegA_now=RegA_trace[i]
                cAMPi_now=cAMPi_trace[i]
                cAMPe_now=cAMPi_trace[i]
                CAR1_now=CAR1_trace[i]
                
                ACA_next,PKA_next,ERK2_next,RegA_next,\
                cAMPi_next,cAMPe_next,CAR1_next=Maeda_agent.update(dt,signal_trace[i])
                
                ACA_trace.append(ACA_next)
                PKA_trace.append(PKA_next)
                ERK2_trace.append(ERK2_next)
                RegA_trace.append(RegA_next)
                cAMPi_trace.append(cAMPi_next)
                # cAMPe_trace.append(cAMPe_next)
                CAR1_trace.append(CAR1_next)
            
            cAMPi_trace = np.array(cAMPi_trace) # make into np array
            cAMPi_trace = cAMPi_trace/Nh_Maeda

            r=np.zeros(NumofCycle-1) # list that stores correlation coefficient to the first peak
            r_new = np.zeros(NumofCycle-1) # correlation coeff take in height difference
            InitPk=cAMPi_trace[int(1*Nt_Maeda/dt): int((1+1*period)*Nt_Maeda/dt)]
            for m in range(NumofCycle-1):
                FollowPk=cAMPi_trace[int((1+(m+1)*period)*Nt_Maeda/dt) : (int((1+(m+1)*period)*Nt_Maeda/dt)+len(InitPk))]
                R = np.corrcoef(InitPk,FollowPk)
                r[m]=R[1,0]
                r_new[m]=r[m]*(1-(sqrt(np.var(InitPk))-sqrt(np.var(FollowPk)))/sqrt(np.var(InitPk)))
#                # check InitPk vs.FollowPk
#                fig3 = plt.figure(figsize=(3,3)); grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
#                ax1= fig3.add_subplot(grid[0, 0])
#                ax1.plot(InitPk); ax1.plot(FollowPk)
#                ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#                plt.show()   
                
            MeanR_Maeda[k,j] = np.mean(r)
            MeanRnew_Maeda[k,j] = np.mean(r_new)
#            # check signal_trace & response trace
#            fig3 = plt.figure(figsize=(4,4)); grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
#            ax1= fig3.add_subplot(grid[0, 0])
#            ax1.plot(t_plot_Maeda,signal_trace)
#            ax1.set_ylabel('extracellular cAMP'); 
#            ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#            ax2= fig3.add_subplot(grid[1:, 0])
#            ax2.plot(t_plot_Maeda,cAMPi_trace)
#            fig3.text(0.5, 0.04, 'E-quality='+str(MeanR_Maeda[k,j])+', new quality=' 
#                      +str(MeanRnew_Maeda[k,j]),fontsize=10, ha='center')
#            plt.show()
            
#%plotperiod- peak width- entrainment quality heatmap
label_font_size = 20
tick_font_size = 16
# fig3 = plt.figure(figsize=(7, 6))
fig3 = plt.figure(figsize=(7, 6))
grid = plt.GridSpec(1, 1, wspace=0.1, hspace=0.1)

ax1= fig3.add_subplot(grid[0, 0])
heatmap = ax1.pcolor(period_space_Maeda, PkWdth_space_Maeda,MeanR_Maeda, cmap='jet')
# heatmap.set_clim(0,1)
fig3.colorbar(heatmap, ax=ax1)
ax1.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Maeda & Loomis 2004 single cellentrainment quality')
plt.show()

#%% Save all outputs in npz file
np.savez('single_cell_entrainment_200218.npz', 
         period_space_Gold=period_space_Gold, PkWdth_space_Gold=PkWdth_space_Gold,
         MeanRnew_Gold=MeanRnew_Gold, MeanR_Gold=MeanR_Gold,
         period_space_Maeda=period_space_Maeda, PkWdth_space_Maeda=PkWdth_space_Maeda,
         MeanRnew_Maeda=MeanRnew_Maeda, MeanR_Maeda=MeanR_Maeda,
         period_space_Sgro=period_space_Sgro, PkWdth_space_Sgro=PkWdth_space_Sgro,
         MeanRnew_Sgro=MeanRnew_Sgro, MeanR_Sgro=MeanR_Sgro,
         period_space_Kamino=period_space_Kamino, PkWdth_space_Kamino=PkWdth_space_Kamino,
         MeanRnew_Kamino=MeanRnew_Kamino, MeanR_Kamino=MeanR_Kamino)
#%% experimental data
Sgro2015Figure4 = scipy.io.loadmat('C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/Figure4Data.mat')
Sgro2015Figure4excel = pd.read_excel(r'C:/Users/ellin/Dropbox/AACP Science/Dicty model review drafts/figures/Sgro2015DataFormattedforPython.xlsx',sheetname='Figure4')
fig3 = plt.figure(figsize=(12,13))
grid = plt.GridSpec(3, 1, wspace=0.2, hspace=0.2)
ax1= fig3.add_subplot(grid[0, 0])

#need to generate the times?
#Fig4Atimes_good = [x*60/271 for x in range(0,271)]
#Fig4Atimes_bad = [x*60/271 for x in range(0,271)]

#ax1.plot(Fig4Atimes_good[0:160],Sgro2015Figure4["goodentrainment"][0,0:160])
#ax1.plot(Fig4Atimes_good[0:160],Sgro2015Figure4["goodentrainment"][1,0:160])
#ax1.plot(Fig4Atimes_good[0:160],Sgro2015Figure4["goodentrainment"][2,0:160])
ax1.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 1 FRET Trace (1 min pulse)"])
ax1.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 2 FRET Trace (1 min pulse)"])
ax1.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 3 FRET Trace (1 min pulse)"])
ax1.set_xlim([0,23])

ax2= fig3.add_subplot(grid[1, 0])
ax2.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 1 FRET Trace (5 min pulse)"])
ax2.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 2 FRET Trace (5 min pulse)"])
ax2.plot(Sgro2015Figure4excel["Time (min)"],Sgro2015Figure4excel["Cell 3 FRET Trace (5 min pulse)"])
ax2.set_xlim([0,23])

#%% load saved npz output file
npzfile = np.load('single_cell_entrainment_200218.npz')
period_space_Gold=npzfile['period_space_Gold']; PkWdth_space_Gold=npzfile['PkWdth_space_Gold']
MeanRnew_Gold= npzfile['MeanRnew_Gold'];  MeanR_Gold=npzfile['MeanR_Gold']
period_space_Maeda=npzfile['period_space_Maeda']; PkWdth_space_Maeda= npzfile['PkWdth_space_Maeda']
MeanRnew_Maeda=npzfile['MeanRnew_Maeda']; MeanR_Maeda=npzfile['MeanR_Maeda']
period_space_Sgro=npzfile['period_space_Sgro']; PkWdth_space_Sgro=npzfile['PkWdth_space_Sgro']
MeanRnew_Sgro= npzfile['MeanRnew_Sgro']; MeanR_Sgro=npzfile['MeanR_Sgro']
period_space_Kamino=npzfile['period_space_Kamino']; PkWdth_space_Kamino=npzfile['PkWdth_space_Kamino']
MeanRnew_Kamino=npzfile['MeanRnew_Kamino']; MeanR_Kamino=npzfile['MeanR_Kamino']

#%% All 4 plots
#% plotperiod- peak width- entrainment quality heatmap

mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']     
#abcd_font_size = 28
#label_font_size=24
#title_font_size = 26
#sublabel_font_size = 22
#trace_width=3
#tick_font_size=20

title_font_size = 22
label_font_size = 20
sublabel_font_size = 20
tick_font_size = 18
trace_width = 3
abcd_font_size = 28


fig3 = plt.figure(figsize=(8,13))
grid = plt.GridSpec(6, 2, wspace=0.6, hspace=1.5)

ax0u= fig3.add_subplot(grid[0, 0])
ax0u.plot(Sgro2015Figure4excel["Time (min)"]-1,Sgro2015Figure4excel["Cell 1 FRET Trace (1 min pulse)"],
                               'k',linewidth=trace_width)
ax0u.plot(Sgro2015Figure4excel["Time (min)"]-1,Sgro2015Figure4excel["Cell 2 FRET Trace (1 min pulse)"],
                               'dimgrey',linewidth=trace_width)
ax0u.plot(Sgro2015Figure4excel["Time (min)"]-1,Sgro2015Figure4excel["Cell 3 FRET Trace (1 min pulse)"],
                               'darkgrey',linewidth=trace_width)
ax0u.set_xlim([0,23]); ax0u.tick_params(axis='both', which='major', labelsize=tick_font_size)
for i in range(3):
    ax0u.axvspan(5+i*6, 5+i*6+1, alpha=0.2, color='b')
ax0u.text(-0.25 , 1.85, 'A', ha='center',va='center',
     transform = ax0u.transAxes, color = 'b', fontsize=abcd_font_size)
ax0u.text(-0.35, -0.5, 'FRET, A.U.',ha='center',va='center', rotation='vertical',
     transform = ax0u.transAxes, color = 'k', fontsize=sublabel_font_size)

#ax0l= fig3.add_subplot(grid[1, 0])
ax0l = fig3.add_axes([0.13, 0.7, 0.305,0.072])
ax0l.plot(Sgro2015Figure4excel["Time (min)"]-1,Sgro2015Figure4excel["Cell 1 FRET Trace (5 min pulse)"],
                               'k',linewidth=trace_width)
ax0l.plot(Sgro2015Figure4excel["Time (min)"]-1,Sgro2015Figure4excel["Cell 2 FRET Trace (5 min pulse)"],
                               'dimgrey',linewidth=trace_width)
ax0l.plot(Sgro2015Figure4excel["Time (min)"]-1,Sgro2015Figure4excel["Cell 3 FRET Trace (5 min pulse)"],
                               'darkgrey',linewidth=trace_width)
ax0l.set_xlim([0,23]); ax0l.tick_params(axis='both', which='major', labelsize=tick_font_size)
for i in range(3):
    ax0l.axvspan(6+i*6, 6+i*6+5, alpha=0.2, color='b')
ax0l.text(0.5, -0.6, 'Time(min)',ha='center',va='center',
     transform = ax0l.transAxes, color = 'k', fontsize=sublabel_font_size)

#ax00= fig3.add_subplot(grid[0:2, 1])
ax00 = fig3.add_axes([0.59, 0.7, 0.31,0.18])
PeriodExp = np.linspace(3, 6, num=4); PkWdthExp = np.linspace(1,5,5)
entrainmentRs = Sgro2015Figure4["entrainmentRs"][:,0:4]; entrainmentRs[entrainmentRs == 0] = 'nan'
heatmap = ax00.pcolor(PeriodExp, PkWdthExp,np.flip(entrainmentRs,0), cmap='jet') 
heatmap.set_clim(0,1)
ax00.set_xlabel('Period(min)', size=sublabel_font_size)
ax00.set_ylabel('Peak Width(min)', size=sublabel_font_size)

cbar=fig3.colorbar(heatmap, ax=ax00,ticks=[0,0.2, 0.4,0.6,0.8,1]);
cbar.ax.tick_params(labelsize = tick_font_size) 
cbar.set_label( 'Entrainment Quality',size=tick_font_size)
ax00.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax00.text(-0.3 , 1.25, 'B', ha='center',va='center',
     transform = ax00.transAxes, color = 'b', fontsize=abcd_font_size)

ax1= fig3.add_subplot(grid[2:4, 0], xticks = [0.8,1.2,1.6], yticks = [0.3,0.6,0.9,1.2])
heatmap = ax1.pcolor(period_space_Gold, PkWdth_space_Gold,MeanR_Gold, cmap='jet') 
# heatmap = ax1.pcolor(period_space_Gold, PkWdth_space_Gold,MeanRnew_Gold, cmap='jet') 
# heatmap.set_clim(0,1)
#cbar=fig3.colorbar(heatmap, ax=ax1, ticks=[0,0.2, 0.4,0.6,0.8,1]);
cbar=fig3.colorbar(heatmap, ax=ax1, ticks=[0.94,0.96,0.98,1]);
cbar.ax.tick_params(labelsize = tick_font_size) 

ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Martiel 1987',color=mycolors[0], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax1.text(-0.3 , 1.1, 'C', ha='center',va='center',
     transform = ax1.transAxes, color = 'g', fontsize=abcd_font_size)

ax2= fig3.add_subplot(grid[2:4, 1], xticks = [0.8,1.2,1.6], yticks = [0.3,0.6,0.9,1.2])
heatmap = ax2.pcolor(period_space_Maeda, PkWdth_space_Maeda,MeanR_Maeda, cmap='jet') # cmap='jet'
heatmap.set_clim(0,1)
cbar=fig3.colorbar(heatmap, ax=ax2, ticks=[0,0.2, 0.4,0.6,0.8,1]);
cbar.ax.tick_params(labelsize = tick_font_size) 

ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Maeda 2004',color=mycolors[1], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax2.text(-0.3 , 1.1, 'D',ha='center',va='center',
     transform = ax2.transAxes, color = 'g', fontsize=abcd_font_size)


ax3= fig3.add_subplot(grid[4:,0], xticks = [0.8,1.2,1.6], yticks = [0.3,0.6,0.9,1.2])
heatmap = ax3.pcolor(period_space_Sgro, PkWdth_space_Sgro,MeanR_Sgro, cmap='jet') # cmap='jet'
heatmap.set_clim(0,1)
cbar=fig3.colorbar(heatmap, ax=ax3, ticks=[0,0.2, 0.4,0.6,0.8,1]);
cbar.ax.tick_params(labelsize = tick_font_size) 

ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Sgro 2015',  color=mycolors[5], fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax3.text(-0.3 , 1.1, 'E', ha='center',va='center',
     transform = ax3.transAxes, color = 'g', fontsize=abcd_font_size)

ax4= fig3.add_subplot(grid[4:,1], xticks = [0.8,1.2,1.6], yticks = [0.3,0.6,0.9,1.2])
heatmap = ax4.pcolor(period_space_Kamino, PkWdth_space_Kamino,MeanR_Kamino, cmap='jet') # cmap='jet'
heatmap.set_clim(0,1)
cbar=fig3.colorbar(heatmap, ax=ax4, ticks=[0,0.2, 0.4,0.6,0.8,1]);
cbar.ax.tick_params(labelsize = tick_font_size) 
#ax4.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
#ax4.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('Kamino 2017', color=mycolors[7],fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
ax4.text(-0.3 , 1.1, 'F', ha='center',va='center',
     transform = ax4.transAxes, color = 'g', fontsize=abcd_font_size)

fig3.text(0.5, 0.04,'Period, A.U.', ha='center', va='center',fontsize=label_font_size)
fig3.text(0.02, 0.4, 'Peak Width, A.U.', ha='center', va='center', rotation='vertical',fontsize=label_font_size)
plt.show()





