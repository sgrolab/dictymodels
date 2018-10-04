% This code simulates spiral waves according to the model by Levine et. al. 
% It calls initial_state.m to set up initial conditions, and Levine1996.m
% to run the simulation. Please make sure to set up the directory to save
% initial conditions (line 11,line 15) and directory to save simulation output videos (line 18)

clear all; close all;

% set up initial conditions
cell_ratio=0.98; % cell density
h=0.2;firing_ratio=0.005;
InitialDir='C:\Users\ellin\Desktop\Levine1996 - Copy\initial cond';% directory to save initial state data
[cell_mask,release_mask,state,X,Y,w,l,x,y]=initial_state(h,cell_ratio,firing_ratio,InitialDir);
%% Start simulation
clear all; clc; 
InitialDir='C:\Users\ellin\Desktop\Levine1996 - Copy\initial cond';% directory that saves initial state data

% set up result output directory
my_dir='C:\Users\ellin\Desktop\Levine1996 - Copy\';

% start simulation
E_min=0.3; % minimum excitability
C_min=4; % minimum exciting threshold
dt_space= 0.228; 
% dt_space=linspace(0.215,0.230,6); % we can scan time step dt
conc_release=62; % rate of cAMP release of fired cells
% conc_release_space=linspace(59,62,6);% we can scan rate of release 
D=0.04; % diffusion coefficient

for m=1:1:length(dt_space)
    dt=dt_space(m);
    t=0:dt:50; % the total time of stimulation 
    FB=[1,2,3]; % excitability feedback, FB=1 (with feedback),2 (without feedback),3 (with uniform feedback)
    % Start simulation
    for n=1:1:length(FB)
        
        [C,E,C_T,state,release_mask]=Levine1996(FB(n),D,E_min,dt,t,conc_release,InitialDir,my_dir);
        disp(strcat('Simulation ended with dt= ',num2str(dt),'D= ',num2str(D),', ConcRel= ',num2str(conc_release),' , FB= ',num2str(FB(n))))
    end
end






