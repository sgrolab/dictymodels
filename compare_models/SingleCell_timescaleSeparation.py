# Single Cell, Timescale Separation 

import sys, pickle 
import numpy as np 
import matplotlib.pyplot as plt 

# add path for function files
sys.path.insert(0,'//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

# import normalization parameters
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/normParams.pickle','rb') as f:
    Nt_Gregor,Nt_Sgro,Nt_Goldbeter,Nt_Maeda,Nt_Kamino,Nh_Gregor,Nh_Sgro,Nh_Goldbeter,Nh_Maeda,Nh_Kamino,Nh_Kamino_offset,Nh_Sgro_offset = pickle.load(f)


#%% Sgro 

import Sgro2015Func as sf

# set parameters 
Params_large={'e':0.1,'tauA':0.09,'tauR':0.9,'g':0.5,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}
Params_small={'e':1,'tauA':0.09,'tauR':0.9,'g':0.5,'c0':1.2,'sigma':0.15,'N':100,\
            'a':0.058,'alpha0':800,'alpha_pde':1000,'Kd':1e-5,'S':1e6,\
            'Nt':27,'Na':3.5,'offset_A':1.5,'flux_thrs':0}

# set time parameters
dt=0.005
t_tot=6*Nt_Sgro
t=np.arange(0,t_tot,dt)

# set initial values 
A0=-1.5
R0=-0.5

# define cAMPe_in trace 
constant_signal= 1
signal_trace=np.zeros(len(t))
signal_trace[int(len(t)/6):] = constant_signal

# set noise amount 
rng = np.random.default_rng(seed=1)
r_large = np.sqrt(dt) * rng.normal(0,1,size=(len(t)))
r_small = np.sqrt(dt) * rng.normal(0,1,size=(len(t)))

# intialize and run cell, large timescacle separation 
cell_large = sf.Cell([1,1],[A0,R0],Params_large,t)
cell_large.run(dt,signal_trace,r_large)

# initialize and run cell, small timescale separation 
cell_small = sf.Cell([0,0],[A0,R0],Params_small,t)
cell_small.run(dt,signal_trace,r_small) 

# calculate nullclines
a = Params_large['a']
Kd = Params_large['Kd']
g = Params_large['g']
c0 = Params_large['c0']

plotAs = np.linspace(-2.5,2.5)
nullclineA_stim = plotAs-1/3*plotAs**3+a*np.log(1+constant_signal/Kd)
nullclineA_noStim = plotAs-1/3*plotAs**3
nullclineR = 1/g*(plotAs+c0)

plt.figure(figsize=(12,8))
plt.subplot(2,2,1)
plt.plot(cell_large.t/Nt_Sgro,cell_large.A)
plt.ylim([-2.5,2.5])
plt.xlim([0,6])
plt.subplot(2,2,2)
plt.plot(cell_small.t/Nt_Sgro,cell_small.A)
plt.ylim([-2.5,2.5])
plt.xlim([0,6])

timeTraceData = [cell_large.t/Nt_Sgro,cell_large.A,cell_small.t/Nt_Sgro,cell_small.A] 

ax = plt.subplot(2,2,3)
plot_cell_large, = ax.plot(cell_large.A,cell_large.R)
plot_cell_small, = ax.plot(cell_small.A,cell_small.R)
plot_nullclineR, = ax.plot(plotAs,nullclineR,color='firebrick',linewidth=2)
plot_nullclineA_noStim, = ax.plot(plotAs,nullclineA_noStim,color='b',linewidth=2)
plot_nullclineA_stim, = ax.plot(plotAs,nullclineA_stim,color='b',linestyle='--',linewidth=2)

xs_cell_large = plot_cell_large.get_ydata()
ys_cell_large = plot_cell_large.get_xdata()
xs_cell_small = plot_cell_small.get_ydata()
ys_cell_small = plot_cell_small.get_xdata()
xs_nullclineR = plot_nullclineR.get_ydata()
ys_nullclineR = plot_nullclineR.get_xdata()
xs_nullclineA_noStim = plot_nullclineA_noStim.get_ydata()
ys_nullclineA_noStim = plot_nullclineA_noStim.get_xdata()
xs_nullclineA_stim = plot_nullclineA_stim.get_ydata()
ys_nullclineA_stim = plot_nullclineA_stim.get_xdata()

nullclineData = [xs_cell_large,ys_cell_large,xs_cell_small,ys_cell_small,xs_nullclineR,ys_nullclineR,xs_nullclineA_noStim,ys_nullclineA_noStim,xs_nullclineA_stim,ys_nullclineA_stim]

# quivers for phase portrait 
A_arr = np.arange(-2.5,3,0.6)
R_arr = np.arange(-1.2,3,0.3)
A_mesh, R_mesh = np.meshgrid(A_arr , R_arr )
dA = A_mesh-(np.power(A_mesh,3))/3-R_mesh+a*np.log(1+constant_signal/Kd)
dR_large = Params_large['e']*(A_mesh-g* R_mesh+c0)
dR_small = Params_small['e']*(A_mesh-g* R_mesh+c0)

phasePortraitData = [A_mesh,R_mesh,dA,dR_large,dR_small]

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_timescaleSeparation_data/timescaleSeparation_Sgro.pickle','wb') as f:
    pickle.dump([timeTraceData,nullclineData,phasePortraitData],f,pickle.HIGHEST_PROTOCOL)

#%% Kamino 

import Kamino2017Func as kf

# set parameters 
tau_large = 1.5 
tau_small = 0.5
n=2
K=4
kt=2
delta=0.01
gamma=3
rho= 0.01
params_large={'tau':tau_large,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}
params_small={'tau':tau_small,'n':n,'K':K,'kt':kt,'delta':delta,\
               'gamma':gamma,'rho':rho}

# set time parameters 
dt=0.001
t_tot=6*Nt_Kamino
t = np.arange(0,t_tot,dt)

# set initial values 
x0=0.01
y0=0.06
z0=0.005

# define cAMPe_in trace 
constant_signal = 1 
signal_trace=np.zeros(len(t))
signal_trace[int(len(t)/6):] = constant_signal

# initialize and run large timescale separation cell 
cell_large = kf.Cell([x0,y0,z0],params_large,t)
cell_large.run(dt,signal_trace)

# initialize and run small timescale separation cell 
cell_small = kf.Cell([x0,y0,z0],params_small,t)
cell_small.run(dt,signal_trace)

timeTraceData = [cell_large.t/Nt_Kamino,cell_large.y,cell_small.t/Nt_Kamino,cell_small.y]

# calculate nullclines 
plotXs = np.linspace(-0.5,2,150)
nullclineX_noStim = 0 + delta
nullclineX_stim = constant_signal + delta
nullclineY_noStim = (0+delta)**n/((0+delta)**n+(np.power(K*plotXs,n)))
nullclineY_stim = (constant_signal+delta)**n/((constant_signal+delta)**n+(np.power(K*plotXs,n)))

nullclineData = [cell_large.x,cell_large.y,cell_small.x,cell_small.y,plotXs,nullclineX_noStim,nullclineX_stim,nullclineY_noStim,nullclineY_stim]

# calculate phase portrait 
x_arr = np.arange(-1,3,0.15)
y_arr = np.arange(-0.5,1.5,0.15)
x_mesh, y_mesh = np.meshgrid(x_arr , y_arr )
dy = (constant_signal+delta)**n/((constant_signal+delta)**n+(np.power(K*x_mesh,n))) - y_mesh
dx_small = 1/tau_small*(constant_signal+delta-x_mesh)
dx_large = 1/tau_large*(constant_signal+delta-x_mesh)
scaleVal = 20

phasePortraitData = [x_mesh,y_mesh,dx_small,dx_large,dy]


with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/sc_timescaleSeparation_data/timescaleSeparation_Kamino.pickle','wb') as f:
    pickle.dump([timeTraceData,nullclineData,phasePortraitData],f,pickle.HIGHEST_PROTOCOL)

