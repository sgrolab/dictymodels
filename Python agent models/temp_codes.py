# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 19:19:49 2019

@author: ellin
"""

# Temp file
#%%
z0First_space_Loomis =np.array([0.1,0.2,0.4,0.8]) #np.array([0.5,10]) # first period extracellular cAMP level # 
FC_space_Loomis=  np.logspace(0.5, 2, num=8) # np.logspace(0.5, 0.5, num=1) 

# Initialize 
PkPrm_Loomis = np.zeros((len(z0First_space_Loomis), len(FC_space_Loomis))) # PkPrm -  peak prominence after the second extracellular cAMP stimulus

k1=1.4; k2=0.9; k3=2.5; k4=1.5; k5=0.6
k6=0.8; k7=2.0; k8=1.3; k9=0.7; k10=1.0
k11=0.3; k12=3.1; k13=1.8; k14=1.5
LaubAgentParam={'k1':k1,'k2':k2,'k3':k3,'k4':k4,'k5':k5,'k6':k6,\
            'k7':k7,'k8':k8,'k9':k9,'k10':k10,'k11':k11,'k12':k12,\
            'k13':k13,'k14':k14}
ACA0=0.1; PKA0=0.1; ERK20=0.1; RegA0=0.1; cAMPi0=0.01; 
cAMPe0=0.1; CAR10=0.1
state0= [ACA0,PKA0,ERK20,RegA0,cAMPi0,cAMPe0,CAR10]

dt=0.0005; t_tot=100; t=list(np.arange(0,t_tot,dt))

#  k_test=np.array([7])
j_test=[3];k_test=[7]

for j in range(len(z0First_space_Loomis)):
    signal_trace=z0First_space_Loomis[j]*np.ones(len(t))
    for k in range(len(FC_space_Gold)):
        
        stim_time_step=int(round(0.5*t_tot/dt)) # at this time second step input is applied
        signal_trace[stim_time_step:] = FC_space_Loomis[k]*z0First_space_Loomis[j]
        
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

        ERK2_trace = np.array(ERK2_trace) # vectorize p_trace
        cAMPi_trace = np.array(cAMPi_trace)
        t_plot_Loomis = np.array(t)
#        # check traces
#        fig3 = plt.figure(figsize=(6, 6))
#        grid = plt.GridSpec(2, 1, wspace=0.3, hspace=0.2)
#        ax1= fig3.add_subplot(grid[0, 0])
#        ax1.plot(t_plot_Loomis,signal_trace)
#        ax1.set_ylabel('extracellular cAMP')
#        ax1.set_title('cAMP from '+str(z0First_space_Loomis[j])+' to FC '+ str(FC_space_Loomis[k]))
#        ax2= fig3.add_subplot(grid[1, 0])
#        ax2.plot(t_plot_Loomis,cAMPi_trace)
#        ax2.set_ylabel(' [cAMP]cyto')
#        ax2.set_ylabel('Time')
#        plt.show()

        end_time_step=len(signal_trace)
        cAMPi_trace_second=cAMPi_trace[stim_time_step:end_time_step] # the second part of trace, second spike
        PkPos, PkProperties = find_peaks(cAMPi_trace_second, prominence=(0,20))
        # Check find_peaks
        # plt.plot(z_trace_later)
        # plt.plot(peaks, z_trace_later[peaks], "x")
        if PkPos.size: # if there is a second spike
            PkPrm_Loomis[j,k]=PkProperties["prominences"][0]
        else:
            PkPrm_Loomis[j,k]=0 # if there is no second spike
                

# PkPrm_Gold_mean=np.mean(PkPrm_Gold,axis=2)
# plot FC vs. second response amplitude
label_font_size = 9
trace_width = 2

fig3 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.2)

ax1= fig3.add_subplot(grid[0, 0])

for i in range(len(z0First_space_Loomis)):
    ax1.plot(FC_space_Loomis,PkPrm_Loomis[i,:], linewidth=trace_width,label='Priming cAMP='+str(z0First_space_Loomis[i]))

ax1.set_ylabel( 'Second spike prominence',fontsize=label_font_size)
ax1.set_xscale('log')
ax1.set_xlabel('Fold change in extracellular cAMP',fontsize=label_font_size)
ax1.set_title('second spike prominence')
leg = ax1.legend();
