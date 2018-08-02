% Sawai 2005
clear all
clc
% initial conditions
h=0.06;firing_ratio=0.005;
[cell_mask,release_mask,state,X,Y,w,l,x,y]=initial_state(h,firing_ratio);
%%
clear all; close all;

% set up directory
global my_dir
my_dir='U:\cilse_research_sgro\Chuqiao\cAMP modeling\Sawai 2005\';
initial_filename='initial_cond\initial_FR_0.005';
load(strcat(my_dir,initial_filename,'.mat'))

global cell_mask X Y w l x y state
global firing_ratio h
h=0.06;firing_ratio=0.005;

% Global constant parameters
global conc_release t_release gamma E_max alpha cell_ratio C_max C_min A D 
conc_release=300; gamma=8; 
E_max=0.93; alpha=0.0005; 
C_max=100; C_min=4; 
cell_ratio=0.98; D=0.00138;
A=(7+2).*(C_max-C_min)./7;

global num_of_dt_release T_ARP T_RRP

% start simulation
E_min=0.5; 
dt_space=0.01; % min
t_release=dt_space; 
beta=0.01;eta=0.0005;

% time step and period, non-dimensional
tic
for m=1:1:length(dt_space)
    dt=dt_space(m);
    t=0:dt:475;
    num_of_dt_release=t_release./dt;
    T_ARP=2./dt; T_RRP=7./dt; % # of time steps in ARP and RRP
    
    % Initializations
    state=repmat(state,1,length(t));
    release_mask=repmat(release_mask,1,length(t));
    C=1.*ones(length(Y),length(t)); % extracellular cAMP concentration
    E=E_min*ones(length(Y),length(t));E(cell_mask==0)=-1; % excitability, mark grids w/o cells as -1
    C_T=C_min*ones(length(Y),length(t)); C_T(cell_mask==0)=-1;% excitation threshold

    % Start simulation
    for n=1:1:1
        [C,E,C_T,state,release_mask]=Sawai2005(eta, beta,dt,t,C,E,C_T,state,release_mask);
    end
end
disp('Simulation ended')
toc;


