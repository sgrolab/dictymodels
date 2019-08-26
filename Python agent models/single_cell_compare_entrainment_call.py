# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:23:10 2019

@author: Chuqiao
"""
# Fold change detection
import numpy as np
import random
import math
import matplotlib.pyplot as plt


from scipy.signal import find_peaks, correlate 

#%% Sgro 2015 with noise
from Sgro2015_agent import Sgro2015_agent
# from TracesXcorr import TracesXcorr

AdptPkPeriod=1.5

#period_space_Sgro =np.linspace(0.9/1.5*AdptPkPeriod,1.7/1.5*AdptPkPeriod,num=7) # period of cAMP stim 
#PkWdth_space_Sgro=  np.linspace(0.3/1.5*AdptPkPeriod, 1.4/1.5*AdptPkPeriod, num=8) # np.logspace(0.5, 0.5, num=1) 
period_space_Sgro =np.linspace(0.5*AdptPkPeriod,1.15*AdptPkPeriod,num=7) # period of cAMP stim 
PkWdth_space_Sgro=  np.linspace(0.2*AdptPkPeriod, 0.95*AdptPkPeriod, num=8) # np.logspace(0.5, 0.5, num=1) 


# Initialize 
MeanR_Sgro = np.empty((len(PkWdth_space_Sgro), len(period_space_Sgro))) 
MeanR_Sgro[:] = np.nan # Mean cross-correlation of spike trace in cycle of stimulus, as nan matrix

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
Nt=27; # normalization factor of t
Na=3.5;  # normalization factor of A
A0=-1.5; R0=-0.5

dt=0.005 
cAMP = 1 # extracellular cAMP
NumofCycle = 8 # cycle stimulation

j_test=[3]; k_test=[5,6]
for j in range(len(period_space_Sgro)): # j_test: #
    period = period_space_Sgro[j]
    for k in  range(len(PkWdth_space_Sgro)): #k_test: #
        PkWdth = PkWdth_space_Sgro[k]
        if PkWdth>=period: 
            pass
        else:
            stim_1cycle=np.zeros(int(period*Nt/dt))
            stim_1cycle[0:int(PkWdth*Nt/dt)]=cAMP;
            signal_trace=np.concatenate((np.zeros(int(1*Nt/dt)), np.tile(stim_1cycle,(NumofCycle)),np.zeros(int(1*Nt/dt))),axis=0);
            
            t = list(np.arange(0,len(signal_trace)*dt,dt))
            t_plot_Sgro=np.array(t)/Nt 
            
            A_trace_orig=[A0]; R_trace_orig=[R0]
            Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)
            for i in range(len(t)-1):
                A_now=A_trace_orig[i]
                R_now=R_trace_orig[i]
                signal_now=signal_trace[i]
        
                A_next,R_next,r_now=Sgro_agent.update(signal_now,dt)
                A_trace_orig.append(A_next)
                R_trace_orig.append(R_next)
       
            A_trace_offset=1.5
            A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
            A_trace_plot=(A_trace_orig+A_trace_offset)/Na;
            
#                        # check signal_trace & response trace
#            fig3 = plt.figure(figsize=(4,4)); grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
#            ax1= fig3.add_subplot(grid[0, 0])
#            ax1.plot(t_plot_Sgro,signal_trace)
#            ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#            ax2= fig3.add_subplot(grid[1:, 0])
#            ax2.plot(t_plot_Sgro,A_trace_plot)
#            plt.show()

            
            r=np.zeros(NumofCycle-1) # list that stores correlation coefficient to the first peak
            InitPk=A_trace_plot[int((1+0*period)*Nt/dt): int((1+1*period)*Nt/dt)]
            # SingleTraces=InitPk
            for m in range(NumofCycle-1):
                FollowPk=A_trace_plot[int((1+(m+1)*period)*Nt/dt) : (int((1+(m+1)*period)*Nt/dt)+len(InitPk))]
                # SingleTraces= np.column_stack((SingleTraces,FollowPk)) 
                for m in range(NumofCycle-1):
                    FollowPk=A_trace_plot[int((1+(m+1)*period)*Nt/dt) : (int((1+(m+1)*period)*Nt/dt)+len(InitPk))]
                    R = np.corrcoef(InitPk,FollowPk)
                    r[m]=R[1,0]
                
#                # check InitPk vs.FollowPk
#                fig3 = plt.figure(figsize=(3,3)); grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
#                ax1= fig3.add_subplot(grid[0, 0])
#                ax1.plot(SecondPk); ax1.plot(FollowPk)
#                ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#                plt.show()    
                    MeanR_Sgro[k,j] = np.mean(r)
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
ax1.set_title('Sgro2015 single cell entrainment quality')
plt.show()

#%% Sgro 2015 withOUT noise
from Sgro2015_agent import Sgro2015_agent
# from TracesXcorr import TracesXcorr
AdptPkPeriod=1.5

period_space_Sgro =np.linspace(0.5*AdptPkPeriod,1.15* AdptPkPeriod,num=7) # period of cAMP stim 
PkWdth_space_Sgro=  np.linspace(0.2*AdptPkPeriod, 0.95*AdptPkPeriod, num=8) # np.logspace(0.5, 0.5, num=1) 

# Initialize 
MeanR_Sgro_NOnoise = np.empty((len(PkWdth_space_Sgro), len(period_space_Sgro))) 
MeanR_Sgro_NOnoise[:] = np.nan # Mean cross-correlation of spike trace in cycle of stimulus, as nan matrix

e=0.1; tauA=0.09; tauR=tauA/e; g=0.5
SgroAgentParam={'e':e,'tauA':tauA,'tauR':tauR,'g':g,'c0':1.2,'sigma':0,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
Nt=27; # normalization factor of t
Na=3.5;  # normalization factor of A
A0=-1.5; R0=-0.5

dt=0.005 
cAMP = 1 # extracellular cAMP
NumofCycle = 8 # cycle stimulation

j_test=[3]; k_test=[5,6]
for j in range(len(period_space_Sgro)): # j_test: #
    period = period_space_Sgro[j]
    for k in  range(len(PkWdth_space_Sgro)): #k_test: #
        PkWdth = PkWdth_space_Sgro[k]
        if PkWdth>=period: 
            pass
        else:
            stim_1cycle=np.zeros(int(period*Nt/dt))
            stim_1cycle[0:int(PkWdth*Nt/dt)]=cAMP;
            signal_trace=np.concatenate((np.zeros(int(1*Nt/dt)), np.tile(stim_1cycle,(NumofCycle)),np.zeros(int(1*Nt/dt))),axis=0);
            
            t = list(np.arange(0,len(signal_trace)*dt,dt))
            t_plot_Sgro=np.array(t)/Nt 
            
            A_trace_orig=[A0]; R_trace_orig=[R0]
            Sgro_agent=Sgro2015_agent([1,1],[A0,R0],SgroAgentParam)
            for i in range(len(t)-1):
                A_now=A_trace_orig[i]
                R_now=R_trace_orig[i]
                signal_now=signal_trace[i]
        
                A_next,R_next,r_now=Sgro_agent.update(signal_now,dt)
                A_trace_orig.append(A_next)
                R_trace_orig.append(R_next)
       
            A_trace_offset=1.5
            A_trace_orig = np.array(A_trace_orig) # vectorize A_trace_orig
            A_trace_plot=(A_trace_orig+A_trace_offset)/Na;
            
#                        # check signal_trace & response trace
#            fig3 = plt.figure(figsize=(4,4)); grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
#            ax1= fig3.add_subplot(grid[0, 0])
#            ax1.plot(t_plot_Sgro,signal_trace)
#            ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#            ax2= fig3.add_subplot(grid[1:, 0])
#            ax2.plot(t_plot_Sgro,A_trace_plot)
#            plt.show()

            
            r=np.zeros(NumofCycle-1) # list that stores correlation coefficient to the first peak
            InitPk=A_trace_plot[int((1+0*period)*Nt/dt): int((1+1*period)*Nt/dt)]
            # SingleTraces=InitPk
            for m in range(NumofCycle-1):
                FollowPk=A_trace_plot[int((1+(m+1)*period)*Nt/dt) : (int((1+(m+1)*period)*Nt/dt)+len(InitPk))]
                # SingleTraces= np.column_stack((SingleTraces,FollowPk)) 
                for m in range(NumofCycle-1):
                    FollowPk=A_trace_plot[int((1+(m+1)*period)*Nt/dt) : (int((1+(m+1)*period)*Nt/dt)+len(InitPk))]
                    R = np.corrcoef(InitPk,FollowPk)
                    r[m]=R[1,0]
                
#                # check InitPk vs.FollowPk
#                fig3 = plt.figure(figsize=(3,3)); grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
#                ax1= fig3.add_subplot(grid[0, 0])
#                ax1.plot(SecondPk); ax1.plot(FollowPk)
#                ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#                plt.show()    
                    MeanR_Sgro_NOnoise[k,j] = np.mean(r)
                    # MeanR_Sgro[k,j] = TracesXcorr(SingleTraces)

#% plotperiod- peak width- entrainment quality heatmap
label_font_size = 20
tick_font_size = 16
fig3 = plt.figure(figsize=(7, 6))
grid = plt.GridSpec(1, 1, wspace=0.1, hspace=0.1)

ax1= fig3.add_subplot(grid[0, 0])
heatmap = ax1.pcolor(period_space_Sgro, PkWdth_space_Sgro,MeanR_Sgro_NOnoise, cmap='jet') # cmap='jet'
heatmap.set_clim(0,1)
fig3.colorbar(heatmap, ax=ax1)

ax1.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Sgro2015 single cell entrainment quality, w/o noise')
plt.show()

#%% Kamino 2017
from Kamino_agent_single_cell import Kamino2017_agent 

AdptPkPeriod=8.5
period_space_Kamino = np.linspace(0.5*AdptPkPeriod,1.15*AdptPkPeriod,num=7) # period of cAMP stim # np.array([5.4,7]) #
PkWdth_space_Kamino=  np.linspace(0.2*AdptPkPeriod,0.95*AdptPkPeriod, num=8) # np.logspace(0.5, 0.5, num=1)  np.array([1.3,4]) # 

# Initialize 
MeanR_Kamino = np.empty((len(PkWdth_space_Kamino), len(period_space_Kamino))) 
MeanR_Kamino[:] = np.nan # Mean cross-correlation of spike trace in cycle of stimulus, as nan matrix

tau=1.5; n=2; K=4; kt=2; delta=0.01
gamma=3; rho= 0.01 # population density, doesn't matter for single cells
KaminoAgentParam={'tau':tau,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
x0=0.01; y0=0.05; z0=0.005

dt=0.001 ; # t_tot=20*Nt; t=list(np.arange(0,t_tot,dt))
cAMP = 1 # extracellular cAMP
NumofCycle = 8 # cycle stimulation

# j_test=[2]; k_test=[3]
for j in range(len(period_space_Kamino)):
    period = period_space_Kamino[j]
    for k in range(len(PkWdth_space_Kamino)):
        PkWdth = PkWdth_space_Kamino[k]
        if PkWdth>=period: 
            pass
        else:
            stim_1cycle=np.zeros(int(period/dt))
            stim_1cycle[0:int(PkWdth/dt)]=cAMP;
            signal_trace=np.concatenate((np.zeros(int(10/dt)), np.tile(stim_1cycle,(NumofCycle)),np.zeros(int(10/dt))),axis=0);
            
            t = list(np.arange(0,len(signal_trace)*dt,dt))
            t_plot_Kamino=np.array(t)
            
            x_trace=[x0]; y_trace=[y0]
            Kamino_agent=Kamino2017_agent([x0,y0,z0],KaminoAgentParam)
            for i in range(len(t)-1):
                x_now=x_trace[i]
                y_now=y_trace[i]
                x_next,y_next,z_next= Kamino_agent.update(1, dt, signal_trace[i])
                x_trace.append(x_next)
                y_trace.append(y_next)
                
       
            y_trace = np.array(y_trace)# vectorize y_trace
            
#            # check signal_trace & response trace
#            fig3 = plt.figure(figsize=(4,4)); grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
#            ax1= fig3.add_subplot(grid[0, 0])
#            ax1.plot(t_plot_Kamino,signal_trace)
#            ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#            ax2= fig3.add_subplot(grid[1:, 0])
#            ax2.plot(t_plot_Kamino,y_trace)
#            plt.show()

            r=np.zeros(NumofCycle-1) # list that stores correlation coefficient to the first peak
            InitPk=y_trace[int(10/dt): int((10+period)/dt)]
            for m in range(NumofCycle-1):
                FollowPk=y_trace[int((10+(m+1)*period)/dt) : (int((10+(m+1)*period)/dt)+len(InitPk))]
                R = np.corrcoef(InitPk,FollowPk)
                r[m]=R[1,0]
#                # check InitPk vs.FollowPk
#                fig3 = plt.figure(figsize=(3,3)); grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
#                ax1= fig3.add_subplot(grid[0, 0])
#                ax1.plot(InitPk); ax1.plot(FollowPk)
#                ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#                plt.show()   
                
            MeanR_Kamino[k,j] = np.mean(r)

#% plotperiod- peak width- entrainment quality heatmap
label_font_size = 20
tick_font_size = 16
# fig3 = plt.figure(figsize=(7, 6))
fig3 = plt.figure(figsize=(7, 6))
grid = plt.GridSpec(1, 1, wspace=0.1, hspace=0.1)

ax1= fig3.add_subplot(grid[0, 0])
heatmap = ax1.pcolor(period_space_Kamino, PkWdth_space_Kamino,MeanR_Kamino, cmap='jet')
# heatmap.set_clim(0,1)
fig3.colorbar(heatmap, ax=ax1)
ax1.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Kamino2017 single cell entrainment quality')
plt.show()


#%% Goldbeter 1987
from Goldbeter1987_agent import Goldbeter1987_agent_3var

AdptPkPeriod = 7.5
period_space_Gold = np.linspace(0.5*AdptPkPeriod,1.15*AdptPkPeriod,num=7) # period of cAMP stim # np.array([5.4,7]) #
PkWdth_space_Gold=  np.linspace(0.2*AdptPkPeriod,0.95*AdptPkPeriod, num=8) # np.logspace(0.5, 0.5, num=1)  np.array([1.3,4]) # 

# Initialize 
MeanR_Gold = np.empty((len(PkWdth_space_Gold), len(period_space_Gold))) 
MeanR_Gold[:] = np.nan # Mean cross-correlation of spike trace in cycle of stimulus, as nan matrix

k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10; L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01; theta=0.01
e= 0.108 # compared to 1
q=4000
sig=0.57 # compared to 0.6
v=12; k= 4 # k prime in the paper
ki=0.958 # compared to 1.7 
kt=0.9
kc=3.58 # compared to 5.4
h=5

Goldbeter3AgentParam={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}

p0=0.8; a0=3; b0=0.9; g0=0

dt=0.001 ; # t_tot=20*Nt; t=list(np.arange(0,t_tot,dt))
cAMP = 1 # extracellular cAMP
NumofCycle = 9 # cycle stimulation

j_test=[3]; k_test=[0,2]

for j in range(len(period_space_Gold)):# j_test:# 
    period = period_space_Gold[j]
    for k in range(len(PkWdth_space_Gold)): # k_test:#
        PkWdth = PkWdth_space_Gold[k]
        if PkWdth>=period: 
            pass
        else:
            stim_1cycle=np.zeros(int(period/dt))
            stim_1cycle[0:int(PkWdth/dt)]=cAMP;
            signal_trace=np.concatenate((np.zeros(int(10/dt)), np.tile(stim_1cycle,(NumofCycle)),np.zeros(int(10/dt))),axis=0);
            
            t = list(np.arange(0,len(signal_trace)*dt,dt))
            t_plot_Gold=np.array(t)
            
            p_trace=[p0]; b_trace=[b0]; g_trace=[g0]
            Goldbeter3_agent=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],Goldbeter3AgentParam)
            for i in range(len(t)-1):
                p_now=p_trace[i]
                p_next,b_next,g_next= Goldbeter3_agent.update(dt,a0,signal_trace[i])
                p_trace.append(p_next)
                b_trace.append(b_next)
                g_trace.append(g_next)
                
       
            b_trace = np.array(b_trace); b_trace = b_trace/np.amax(b_trace)
            
#            # check signal_trace & response trace
#            fig3 = plt.figure(figsize=(4,4)); grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
#            ax1= fig3.add_subplot(grid[0, 0])
#            ax1.plot(t_plot_Gold,signal_trace)
#            ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#            ax2= fig3.add_subplot(grid[1:, 0])
#            ax2.plot(t_plot_Gold,b_trace)
#            plt.show()

            r=np.zeros(NumofCycle-1) # list that stores correlation coefficient to the first peak
            InitPk=b_trace[int(10/dt): int((10+period)/dt)]
            for m in range(NumofCycle-1):
                FollowPk=b_trace[int((10+(m+1)*period)/dt) : (int((10+(m+1)*period)/dt)+len(InitPk))]
                R = np.corrcoef(InitPk,FollowPk)
                r[m]=R[1,0]
#                # check InitPk vs.FollowPk
#                fig3 = plt.figure(figsize=(3,3)); grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
#                ax1= fig3.add_subplot(grid[0, 0])
#                ax1.plot(InitPk); ax1.plot(FollowPk)
#                ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#                plt.show()   
                
            MeanR_Gold[k,j] = np.mean(r)

#% plotperiod- peak width- entrainment quality heatmap
label_font_size = 20
tick_font_size = 16
# fig3 = plt.figure(figsize=(7, 6))
fig3 = plt.figure(figsize=(7, 6))
grid = plt.GridSpec(1, 1, wspace=0.1, hspace=0.1)

ax1= fig3.add_subplot(grid[0, 0])
heatmap = ax1.pcolor(period_space_Gold, PkWdth_space_Gold,MeanR_Gold, cmap='jet')
# heatmap.set_clim(0.8,1)
fig3.colorbar(heatmap, ax=ax1)
ax1.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Goldbeter1986 single cell entrainment quality')
plt.show()

        
   
#%% Laub and Loomis 1998
from LaubLoomis1998_agent import LaubLoomis1998_agent

AdptPkPeriod=6
period_space_Laub = np.linspace(0.5*AdptPkPeriod,1.15*AdptPkPeriod,num=7) # period of cAMP stim # np.array([5.4,7]) #
PkWdth_space_Laub=  np.linspace(0.2*AdptPkPeriod,0.95*AdptPkPeriod, num=8) # np.logspace(0.5, 0.5, num=1)  np.array([1.3,4]) # 

# Initialize 
MeanR_Laub = np.empty((len(PkWdth_space_Laub), len(period_space_Laub))) 
MeanR_Laub[:] = np.nan # Mean cross-correlation of spike trace in cycle of stimulus, as nan matrix

k1=1.4; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=2.0; k8=1.3; k9=0.7; k10=1.0
k11=0.3; k12=3.1; k13=1.8; k14=1.5
LaubAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}
ACA0=0.1; PKA0=0.1; ERK20=0.1; RegA0=0.1; cAMPi0=0.01; 
cAMPe0=0.1; CAR10=0.1
state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]

dt=0.0005
cAMP = 1 # extracellular cAMP
NumofCycle = 4 # cycle stimulation

j_test=[3]; k_test=[0,2]

for j in range(len(period_space_Gold)):# j_test:# 
    period = period_space_Laub[j]
    for k in range(len(PkWdth_space_Laub)): # k_test:#
        PkWdth = PkWdth_space_Laub[k]
        if PkWdth>=period: 
            pass
        else:
            stim_1cycle=np.zeros(int(period/dt))
            stim_1cycle[0:int(PkWdth/dt)]=cAMP;
            signal_trace=np.concatenate((np.zeros(int(10/dt)), np.tile(stim_1cycle,(NumofCycle)),np.zeros(int(10/dt))),axis=0);
            
            t = list(np.arange(0,len(signal_trace)*dt,dt))
            t_plot_Laub=np.array(t)
            
            ACA_trace=[ACA0]; PKA_trace=[PKA0]; ERK2_trace=[ERK20]; RegA_trace=[RegA0]; 
            cAMPi_trace=[cAMPi0]; cAMPe_trace=[cAMPe0]; CAR1_trace=[CAR10]
            LaubLoomis_agent=LaubLoomis1998_agent([1,1],state0,LaubAgentParam)
            for i in range(len(t)-1):
                ACA_now=ACA_trace[i]
                PKA_now=PKA_trace[i]
                ERK2_now=ERK2_trace[i]
                RegA_now=RegA_trace[i]
                cAMPi_now=cAMPi_trace[i]
                cAMPe_now=cAMPi_trace[i]
                CAR1_now=CAR1_trace[i]
                
                ACA_next,PKA_next,ERK2_next,RegA_next,\
                cAMPi_next,cAMPe_next,CAR1_next=LaubLoomis_agent.update(dt,signal_trace[i])
                
                ACA_trace.append(ACA_next)
                PKA_trace.append(PKA_next)
                ERK2_trace.append(ERK2_next)
                RegA_trace.append(RegA_next)
                cAMPi_trace.append(cAMPi_next)
                # cAMPe_trace.append(cAMPe_next)
                CAR1_trace.append(CAR1_next)
                
       
            ERK2_trace = np.array(ERK2_trace) # vectorize traces
            cAMPi_trace = np.array(cAMPi_trace)
            
#            # check signal_trace & response trace
#            fig3 = plt.figure(figsize=(4,4)); grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)
#            ax1= fig3.add_subplot(grid[0, 0])
#            ax1.plot(t_plot_Laub,signal_trace)
#            ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#            ax2= fig3.add_subplot(grid[1:, 0])
#            ax2.plot(t_plot_Laub,cAMPi_trace)
#            plt.show()

            r=np.zeros(NumofCycle-1) # list that stores correlation coefficient to the first peak
            InitPk=cAMPi_trace[int(10/dt): int((10+period)/dt)]
            for m in range(NumofCycle-1):
                FollowPk=cAMPi_trace[int((10+(m+1)*period)/dt) : (int((10+(m+1)*period)/dt)+len(InitPk))]
                R = np.corrcoef(InitPk,FollowPk)
                r[m]=R[1,0]
                
#                # check InitPk vs.FollowPk
#                fig3 = plt.figure(figsize=(3,3)); grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)
#                ax1= fig3.add_subplot(grid[0, 0])
#                ax1.plot(InitPk); ax1.plot(FollowPk)
#                ax1.set_ylabel('extracellular cAMP'); ax1.set_title('period '+str(period)+' stim peak width '+ str(PkWdth))
#                plt.show()   
                
            MeanR_Laub[k,j] = np.mean(r)

#% plotperiod- peak width- entrainment quality heatmap
label_font_size = 20
tick_font_size = 16
# fig3 = plt.figure(figsize=(7, 6))
fig3 = plt.figure(figsize=(7, 6))
grid = plt.GridSpec(1, 1, wspace=0.1, hspace=0.1)

ax1= fig3.add_subplot(grid[0, 0])
heatmap = ax1.pcolor(period_space_Laub, PkWdth_space_Laub,MeanR_Laub, cmap='jet')
# heatmap.set_clim(0,1)
fig3.colorbar(heatmap, ax=ax1)
ax1.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
ax1.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Laub & Loomis 1998 single cell entrainment quality')
plt.show()

#%% All 4 plots
#% plotperiod- peak width- entrainment quality heatmap
title_font_size = 20 
label_font_size = 20
tick_font_size = 16

fig3 = plt.figure(figsize=(12,13))
grid = plt.GridSpec(2, 2, wspace=0.2, hspace=0.3)

ax1= fig3.add_subplot(grid[0, 0])
# heatmap = ax1.pcolor(period_space_Gold*10/np.amax(period_space_Gold), PkWdth_space_Gold* 1.4/1.7*10/np.amax(PkWdth_space_Gold),MeanR_Gold, cmap='jet') # cmap='jet'
heatmap = ax1.pcolor(period_space_Gold, PkWdth_space_Gold,MeanR_Gold, cmap='jet') 
# heatmap.set_clim(0,1)
cbar=fig3.colorbar(heatmap, ax=ax1);cbar.ax.tick_params(labelsize = tick_font_size) 
#ax1.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
#ax1.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1.set_title('Goldbeter 1987', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax2= fig3.add_subplot(grid[0, 1])
heatmap = ax2.pcolor(period_space_Laub, PkWdth_space_Laub,MeanR_Laub, cmap='jet') # cmap='jet'
# heatmap.set_clim(0,1)
cbar=fig3.colorbar(heatmap, ax=ax2); cbar.ax.tick_params(labelsize = tick_font_size) 
#ax2.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
#ax2.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.set_title('Laub & Loomis 1998', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax3= fig3.add_subplot(grid[1,0])
heatmap = ax3.pcolor(period_space_Sgro, PkWdth_space_Sgro,MeanR_Sgro, cmap='jet') # cmap='jet'
heatmap.set_clim(0,1)
cbar=fig3.colorbar(heatmap, ax=ax3);cbar.ax.tick_params(labelsize = tick_font_size) 
#ax3.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
#ax3.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax3.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax3.set_title('Sgro 2015', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

ax4= fig3.add_subplot(grid[1,1])
heatmap = ax4.pcolor(period_space_Kamino, PkWdth_space_Kamino,MeanR_Kamino, cmap='jet') # cmap='jet'
# heatmap.set_clim(0,1)
cbar=fig3.colorbar(heatmap, ax=ax4);cbar.ax.tick_params(labelsize = tick_font_size) 
#ax4.set_ylabel('Peak Width, A.U.',fontsize=label_font_size)
#ax4.set_xlabel('Entrainment period, A.U.',fontsize=label_font_size)
ax4.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax4.set_title('Kamino & Sawai 2017', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

fig3.text(0.5, 0.04,'Entrainment period, A.U.', ha='center', va='center',fontsize=label_font_size)
fig3.text(0.06, 0.5, 'Peak Width, A.U.', ha='center', va='center', rotation='vertical',fontsize=label_font_size)

plt.show()





