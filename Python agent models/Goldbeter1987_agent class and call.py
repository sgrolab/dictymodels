# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 16:17:04 2019

@author: Chuqiao Huyan

Goldbeter 1987 model
"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt

#%% Define 4 variable class
class Goldbeter1987_agent_4var:
    def __init__(self,pos,state,param):
        self.pos=pos
        self.state=state
        self.param=param
        self.p_now=state[0]
        self.a_now=state[1]
        self.b_now=state[2] 
        self.g_now=state[3] # initial state input as a list variable [p0, b0, g0]
        
    
    def update(self,dt,signals='none'):
        k1=self.param['k1']   
        k2=self.param['k2']  
        L1=self.param['L1']  
        L2=self.param['L2']  
        c=self.param['c']  
        lamda=self.param['lamda']  
        theta=self.param['theta']  
        e=self.param['e']  
        q=self.param['q'] 
        sig=self.param['sig']  
        v=self.param['v']  
        k=self.param['k']  
        ki=self.param['ki']  
        kt=self.param['kt']  
        kc=self.param['kc']  
        h=self.param['h']  
        
        f1 = (k1+k2*self.g_now)/(1+self.g_now)
        f2 = (k1*L1+k2*L2*c*self.g_now)/(1+c*self.g_now)
        Y= self.p_now*self.g_now/(1+ self.g_now)
        Ysq = (Y)**2
        PI = self.a_now*(lamda*theta + e*Ysq)/(1 + theta*self.a_now + (1 + self.a_now)*e*Ysq)
        
        # p_next=0.8
        p_next = self.p_now+dt*(-f1*self.p_now+f2*(1-self.p_now))
        a_next = self.a_now+dt*(v-k*self.a_now-sig*PI)
        b_next = self.b_now+dt*(q*sig*PI-(ki+kt)*self.b_now)
        
        if isinstance(signals, str):
            g_next = self.g_now+dt*(kt/h*self.b_now-kc*self.g_now)
        else: 
            g_next = signals # if there is contant cAMP input to the system
        
        self.p_now = p_next
        self.a_now = a_next
        self.b_now = b_next
        self.g_now = g_next
        return p_next,a_next,b_next,g_next
    
#%% parameters
        
k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10
L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01
theta=0.01
e=1
q=4000
sig=0.6
v=12
k= 4        # k prime in the paper
ki=1.7 
kt=0.9
kc=5.4
h=5

param={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}

p0=0.6
a0=3
b0=0.9
g0=0

# Fig 2, 4 variable model, autonomous oscillation
test_agent=Goldbeter1987_agent_4var([1,1],[p0,a0,b0,g0],param)
# initializations
p_trace=[p0]
a_trace=[a0]
b_trace=[b0]
g_trace=[g0]


dt=0.001
t_tot=35
t=list(np.arange(0,t_tot,dt))

for i in range(len(t)-1):
    p_now=p_trace[i]
    
    p_next,a_next,b_next,g_next, Y_now= test_agent.update(dt)
    p_trace.append(p_next)
    a_trace.append(a_next)
    b_trace.append(b_next)
    g_trace.append(g_next)
    Y_trace.append(Y_now)
   
# Plot the results
p_trace = np.array(p_trace) # vectorize p_trace
deltaT_trace = 1 - p_trace
a_trace = np.array(a_trace)
b_trace = np.array(b_trace)
g_trace = np.array(g_trace)
Y_trace = p_trace*g_trace/(1+g_trace)+(1-p_trace)*param['c']*g_trace/(1+param['c']*g_trace)

t_plot = np.array(t)

fig,ax=plt.subplots(2,1)

ax[0].plot(t_plot,p_trace,label='pT')
ax[0].plot(t_plot,a_trace, label=r'$\alpha$') 
ax[0].plot(t_plot,b_trace/60, label=r'$\beta$') # scale to plot on same axis
ax[0].plot(t_plot,g_trace/60, label=r'$\gamma$') # scale to plot on same axis
ax[0].set_title('Fig 2a')
ax[0].set_xlabel('Time(min)')
ax[0].set_ylabel('Levels')
leg0 = ax[0].legend()

ax[1].plot(t_plot, Y_trace, label='Y, receptor saturation' )
ax[1].plot(t_plot, g_trace, label=r'$\gamma$')
ax[1].plot(t_plot,deltaT_trace,label=r'$\delta$T')
ax[1].set_title('Fig 2b')
ax[1].set_xlabel('Time(min)')
ax[1].set_ylabel('Levels')
leg1 = ax[1].legend()
plt.show()

#%% Define 3 variable class agent

class Goldbeter1987_agent_3var(Goldbeter1987_agent_4var):
    def update(self,dt, a, signals='none'):
        k1=self.param['k1']   
        k2=self.param['k2']  
        L1=self.param['L1']  
        L2=self.param['L2']  
        c=self.param['c']  
        lamda=self.param['lamda']  
        theta=self.param['theta']  
        e=self.param['e']  
        q=self.param['q'] 
        sig=self.param['sig']  
        v=self.param['v']  
        k=self.param['k']  
        ki=self.param['ki']  
        kt=self.param['kt']  
        kc=self.param['kc']  
        h=self.param['h']  
        
        f1 = (k1+k2*self.g_now)/(1+self.g_now)
        f2 = (k1*L1+k2*L2*c*self.g_now)/(1+c*self.g_now)
        Ysq = (self.p_now*self.g_now/(1+ self.g_now))**2
        PI = a*(lamda*theta + e*Ysq)/(1 + theta*a + (1 + a)*e*Ysq)
        
        # p_next=0.8
        p_next = self.p_now+dt*(-f1*self.p_now+f2*(1-self.p_now))
        b_next = self.b_now+dt*(q*sig*PI-(ki+kt)*self.b_now)
        
        if isinstance(signals, str):
            g_next = self.g_now+dt*(kt/h*self.b_now-kc*self.g_now)
        else: 
            g_next = signals # if there is contant cAMP input to the system
        
        self.p_now = p_next
        self.b_now = b_next
        self.g_now = g_next
        
        return p_next, b_next, g_next
#%% parameters
        
k1 = 0.036     # per min
k2 = 0.666    # per min
L1 = 10
L2 = 0.005 
c = 10;           # 0.15 ~ 50
lamda=0.01
theta=0.01
e=0.108 # compared to 1
q=4000
sig=0.57 # compared to 0.6
v=12
k= 4        # k prime in the paper
ki=0.958 # compared to 1.7 
kt=0.9
kc=3.58 # compared to 5.4
h=5

param={'k1':k1,'k2':k2,'L1':L1,'L2':L2, 'c':c, 'lamda':lamda,\
            'theta':theta, 'e':e, 'q':q,'sig':sig, 'v':v, 'k':k, \
            'ki':ki,'kt':kt, 'kc':kc,'h':h}

p0=0.8
a0=3
b0=0.9
g0=0

# Fig 2, 4 variable model, autonomous oscillation
agent_3var=Goldbeter1987_agent_3var([1,1],[p0,a0,b0,g0],param)
# initializations
p_trace=[p0]
b_trace=[b0]
g_trace=[g0]


dt=0.001
t_tot=25
t=list(np.arange(0,t_tot,dt))

stim_time_step=3000 # at this time step step input is applied
constant_signal=0.3;
signal_trace=np.zeros(len(t))
signal_trace[stim_time_step:] = constant_signal


for i in range(len(t)-1):
    p_now=p_trace[i]
    if i< 3000:
        p_next,b_next,g_next= agent_3var.update(dt,a0,'none')
    else:
        p_next,b_next,g_next= agent_3var.update(dt,a0,constant_signal)
    p_trace.append(p_next)
    b_trace.append(b_next)
    g_trace.append(g_next)
        
   
# Convert into np array
p_trace = np.array(p_trace) # vectorize p_trace
b_trace = np.array(b_trace)

t_plot = np.array(t)

#%% Fig 5
label_font_size = 25
trace_width = 6.0

fig5 = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(3, 1, wspace=0.3, hspace=0.2)

ax1= fig5.add_subplot(grid[0, 0])
ax1.plot(t_plot,signal_trace, linewidth=trace_width)
ax1.set_ylabel( r'$cAMP_{ext}$ input' ,fontsize=label_font_size)
# ax1.set_xlabel('Time',fontsize=label_font_size)

ax2= fig5.add_subplot(grid[1:, 0])
line1=ax2.plot(t_plot,b_trace, color='g',linewidth=trace_width)
ax2.set_ylabel(r'$cAMP_{cyto},  \beta$',fontsize=label_font_size)
ax2.yaxis.label.set_color('g')

ax3 = fig5.add_subplot(grid[1:, 0], sharex=ax2, frameon=False)
line2=ax3.plot(t_plot,p_trace, color='b', linewidth=trace_width) # right axis
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")
ax3.set_ylabel( r'Active form receptor, $\rho_T$' ,fontsize=label_font_size)
ax3.yaxis.label.set_color('b')
# ax3.legend((line1, line2), (r'$\beta$',r'$\rho_T$'))
ax3.set_xlabel('Time',fontsize=label_font_size)

plt.show()
