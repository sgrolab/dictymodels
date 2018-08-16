% Sawai 2005
clear all
clc

% setting output directory
global my_dir
my_dir='U:\cilse_research_sgro\Chuqiao\cAMP modeling\Sawai 2005\';

% Generating initial conditions and save it in my_dir
h=0.06;firing_ratio=0.005;
[cell_mask,release_mask,state,X,Y,w,l,x,y,select]=initial_state(my_dir,h,firing_ratio);
%% 
clear all

% setting output directory
global my_dir
my_dir='U:\cilse_research_sgro\Chuqiao\cAMP modeling\Sawai 2005\';


global cell_mask X Y w l x y 
initial_filename='initial_cond\initial_FR_0.005xy75';
load(strcat(my_dir,initial_filename,'.mat'))
release_mask_prime=release_mask;release_mask_prime(select)=1;

global firing_ratio h
h=0.06;firing_ratio=0.005;

% Global constant parameters
global conc_release t_release gamma E_max  cell_ratio C_max C_min A D E_min C0
conc_release=200; gamma=8; % /min
E_max=0.93; 
C_max=100; C_min=4; 
cell_ratio=0.98; D=0.02;%0.0138; % mm2/min
A=(7+2).*(C_max-C_min)./7;
C0=0.1;
global num_of_dt_release T_ARP T_RRP

% start simulation
E_min=0.5; 
dt=0.01; % min
t_release=dt; % release happens in 1 time step
beta=0.01;eta=0.0005;

% initial conditions 
E0=E_min*ones(length(Y),1);E0(cell_mask==0)=-1;
C_T0=C_min*ones(length(Y),1);C_T0(select)=C_max;C_T0(cell_mask==0)=-1;
release_mask0=release_mask;
state0=state;

% start simulation 
tic

t=0:dt:25;
num_of_dt_release=t_release./dt; % cAMP released in 1 time step
T_ARP=2./dt; T_RRP=7./dt; % # of time steps in ARP and RRP

% Start simulation
for n=1:1:1
    % Initializations
    C=zeros(length(Y),length(t)); C(:,1)=C0.*ones(length(Y),1); % extracellular cAMP concentration
    E=zeros(length(Y),length(t));E(:,1)=E0; % excitability, mark grids w/o cells as -1
    C_T=zeros(length(Y),length(t));C_T(:,1)=C_T0; % excitation threshold
    release_mask=zeros(length(Y),length(t));release_mask(:,1)=release_mask0;
    state=zeros(length(Y),length(t));state(:,1)=state0;
    [C,E,C_T,state,release_mask]=Sawai2005(eta, beta,dt,t,C,E,C_T,state,release_mask);
end

disp('Simulation ended')
toc;


