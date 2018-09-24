% Levine 1996
% initial conditions
h=0.2;firing_ratio=0.005;
[cell_mask,release_mask,state,X,Y,w,l,x,y]=initial_state(h,firing_ratio);
%%
clear all; close all;
global cell_mask X Y w l x y state
global firing_ratio h
h=1;firing_ratio=0.005;
% set up directory
global my_dir
my_dir='U:\cilse_research_sgro\Chuqiao\cAMP modeling\Levine 1996\';
initial_filename='initial_FR_0.005_1';
load(strcat(my_dir,'initial_cond\',initial_filename,'.mat'))

%% Global constant parameters
global conc_release t_release gamma E_max alpha beta cell_ratio C_max C_min A D
conc_release=300; t_release=1; gamma=8; 
E_max=0.93; alpha=0.0005; beta=1.24;
C_max=100; C_min=4; 
cell_ratio=0.98; D=1;

% uniform feedback dE each time step
global dE
global num_of_dt_release T_ARP T_RRP

% start simulation
E_min=0.3; 
dt_space=0.22;% 0.18:0.01:0.25; % scan parameter space

% time step and period, non-dimensional

for m=1:1:length(dt_space)
    dt=dt_space(m);
    t=0:dt:100;
    num_of_dt_release=t_release./dt;
    T_ARP=8./dt; T_RRP=2./dt; % # of time steps in ARP and RRP
    % T_ARP=2./dt; T_RRP=7./dt; % Sawai 2015
    A=(T_RRP+T_ARP).*(C_max-C_min)./T_RRP;
    % uniform feedback 
    dE=E_max-E_min/(200./dt);
    
    % Initializations
    C=1.*ones(length(Y),1); % extracellular cAMP concentration
    % release_mask(cell_mask==1)=1;
    E=E_min*ones(length(Y),1);E(cell_mask==0)=-1; % excitability, mark grids w/o cells as -1
    C_T=C_min*ones(length(Y),1); C_T(cell_mask==0)=-1;% excitation threshold

%     figure % show the cell initial states
%     surf(x,y,reshape(state,[w,l]))
%     shading interp; view(2);
    FB=[2,3]; % feedback=1 (with feedback),2 (without feedback),3 (with uniform feedback)
    % Start simulation
    for n=1:1: length(FB)
        [C,E,C_T,state,release_mask]=Levine1996(FB(n),E_min,dt,t,C,E,C_T,state,release_mask);
    end
end






