% Sgro 2015 single cell
clear all;clc;
% parameters
param_WT.e=0.1; %tauA/tauR; %episilon
param_WT.tauA=0.09;
param_WT.tauR=param_WT.tauA/param_WT.e;

param_WT.g=0.5; %gamma
param_WT.c0=1.2;
param_WT.sigma=0.15; %sigma, noise strength

param_WT.N=100; % number of cells in population simulations
param_WT.a=0.058;
param_WT.alpha0=800;
param_WT.alpha_pde= 1000;
param_WT.Kd= 1e-5;
param_WT.S=1e6;                          

param_WT.Nt=27; % normalization factor of t
param_WT.Na=3.5;  % normalization factor of A
param_WT.offset_A=1.5;

% PKA mutant
param_PKAR.e=0.1; %tauA/tauR; %episilon
param_PKAR.tauA=0.09;
param_PKAR.tauR=param_PKAR.tauA/param_PKAR.e;

a1=0.61; b1=1.25;
param_PKAR.g=1/a1; 
param_PKAR.c0=b1/a1;
param_PKAR.sigma=0.15; %sigma, noise strength

param_PKAR.N=100; % number of cells in population simulations
param_PKAR.a=0.058;
param_PKAR.alpha0=800;
param_PKAR.alpha_pde= 1000;
param_PKAR.Kd= 1e-5;
param_PKAR.S=1e6;                          

param_PKAR.Nt=27; % normalization factor of t
param_PKAR.Na=3.5;  % normalization factor of A
param_PKAR.offset_A=1.5;

%% Step stimulus, fig 1 G& H
dt=0.005;t=0:dt:15*param_WT.Nt;
t_plot=t./param_WT.Nt;
cAMP=1;
StimInitTime=2.5; % stimulation time after Nt normalization 
stim=zeros(size(t));stim(StimInitTime/dt*param_WT.Nt:end)=cAMP;
A0=-1.5;R0=-0.5;

[tout,xout]=euler_solver(@Sgro_single_cell_new,t,[A0;R0],param_WT,stim);
A_WT_orig=xout(1,:);R_WT_orig=xout(2,:); % un-scaled, un-offset A and R after simulation
A_WT=(A_WT_orig+param_WT.offset_A)./param_WT.Na;
R_WT=R_WT_orig./param_WT.Na;

% Nullclines
A_null_WT=linspace(-2.5,2.5,200);
dAdt_null_no_stim_WT=A_null_WT-1/3*A_null_WT.^3;% dA/dt=0 nullcline w/o stim
dAdt_null_stim_WT=A_null_WT-1/3*A_null_WT.^3+param_WT.a*log(1+cAMP/param_WT.Kd)*ones(size(A_null_WT));% dA/dt=0 nullcline with stim
dRdt_null_WT=1/param_WT.g*(A_null_WT+param_WT.c0*ones(size(A_null_WT))); % dR/dt=0 nullcline

% traces and nullclines for PKAR mutant
[tout,xout]=euler_solver(@Sgro_single_cell_new,t,[A0;R0],param_PKAR,stim);
A_PKAR_orig=xout(1,:);R_PKAR_orig=xout(2,:); % un-scaled, un-offset A and R after simulation
A_PKAR=(A_PKAR_orig+param_PKAR.offset_A)./param_PKAR.Na;
R_PKAR=R_PKAR_orig./param_PKAR.Na;

A_null_PKAR=linspace(-2.5,2.5,200);
dAdt_null_no_stim_PKAR=A_null_PKAR-1/3*A_null_PKAR.^3;% dA/dt=0 nullcline w/o stim
dAdt_null_stim_PKAR=A_null_PKAR-1/3*A_null_PKAR.^3+param_PKAR.a*log(1+cAMP/param_PKAR.Kd)*ones(size(A_null_PKAR));% dA/dt=0 nullcline with stim
dRdt_null_PKAR=1/param_PKAR.g*(A_null_PKAR+param_PKAR.c0*ones(size(A_null_PKAR))); % dR/dt=0 nullcline

figure();  
subplot(2,2,1)
plot(A_null_WT,dRdt_null_WT,'r','LineWidth',1)
hold on
plot(A_null_WT,dAdt_null_no_stim_WT,'Color',[0 0.5 0],'LineWidth',1)
hold on 
plot(A_null_WT,dAdt_null_stim_WT,'k','LineWidth',1)
ylim([-2 2.5])
hold on
scatter(A_WT_orig(1:200:end),R_WT_orig(1:200:end),30,'.','k')

subplot(2,2,2) 
plot(t_plot,A_WT,'color',[0.1 0.5 0.25],'LineWidth',1.5); hold on
% line([StimInitLine StimInitLine],'Color',[0 1 0],'LineWidth',1.5,'LineStyle','--')
ylabel('Amplitude'); xlabel('Time (T)')
title(['WT,','cAMP=',num2str(cAMP),', \gamma=',num2str(param_WT.g),', c0=',num2str(param_WT.c0),', \epsilon=',num2str(param_WT.e),])
set(gca,'FontSize',10)

% PKAR mutant
subplot(2,2,3)
plot(A_null_PKAR,dRdt_null_PKAR,'r','LineWidth',1)
hold on
plot(A_null_PKAR,dAdt_null_no_stim_PKAR,'Color',[0 0.5 0],'LineWidth',1)
hold on 
plot(A_null_PKAR,dAdt_null_stim_PKAR,'k','LineWidth',1)
ylim([-2 2.5])
hold on
scatter(A_PKAR_orig(1:200:end),R_PKAR_orig(1:200:end),30,'.','k')

subplot(2,2,4) 
plot(t_plot,A_PKAR,'color',[0.1 0.5 0.25],'LineWidth',1.5); hold on
% line([StimInitLine StimInitLine],'Color',[0 1 0],'LineWidth',1.5,'LineStyle','--')
ylabel('Amplitude'); xlabel('Time (T)')
title(['PKAR-,','cAMP=',num2str(cAMP),', a1=',num2str(a1),', b1=',num2str(b1),', \epsilon=',num2str(param_PKAR.e),])
set(gca,'FontSize',10)

OutFolderName='E:\bu\research\cAMP modeling\PKA mutant modeling\20190125_gamma_c0_single cell\';
saveas(gcf,[OutFolderName,'a1_',num2str(a1),'_b1_',num2str(b1),'_e_',num2str(param_WT.e),'_cAMP_',num2str(cAMP),'.png'])
close all

%% Fig 2C- initial peak width
logcAMP=0:0.2:4.4;
n=numel(logcAMP);
initPKWdth_PKAR=zeros(1,n);
initPKWdth_WT=zeros(1,n);
dt=0.005; t=0:dt:200;
A0=-1.5;R0=-0.5;
StimInitTime=1;
for i=1:n % n=1 to test findpeaks 
    cAMP=10.^(logcAMP(i));
    stim=zeros(size(t));stim(StimInitTime/dt*param_WT.Nt:end)=cAMP;
    % PKAR mutant
    [tout,xout_PKAR]=euler_solver(@Sgro_single_cell_new,t,[A0;R0],param_PKAR,stim);
    A_PKAR=xout_PKAR(1,:);
    
    yy_PKAR = smooth(A_PKAR,2000,'moving'); % moving average filter
    
%     plot(t,A_PKAR,'color',[0 1 0],'LineWidth',3); hold on % test findpeaks
%     findpeaks(yy_PKAR,t,'Annotate','extents','WidthReference','halfheight',...
%         'MinPeakHeight',0.5,'MinPeakDistance',1,'MinPeakProminence',1);hold on 
    
    [PKS_PKAR,LOCS_PKAR,W_PKAR]=findpeaks(yy_PKAR,t,'WidthReference','halfheight',...
        'MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',1);
    if isempty(W_PKAR)
        initPKWdth_PKAR(i)=NaN;
    else
        initPKWdth_PKAR(i)=W_PKAR(1);
    end
    % WT
    i=23;
        [tout,xout_WT]=euler_solver(@Sgro_single_cell_new,t,[A0;R0],param_WT,stim);
    A_WT=xout_WT(1,:);
    yy_WT= smooth(A_WT,200,'moving'); % moving average filter
%     
%     plot(t,A_WT,'color',[0 1 0],'LineWidth',3); hold on % test findpeaks
%     findpeaks(yy_WT,t,'Annotate','extents','WidthReference','halfheight','MinPeakHeight',0.5,...
%         'MinPeakDistance',0.5,'MinPeakProminence',0.5);hold on
    
    [PKS_WT,LOCS_WT,W_WT]=findpeaks(yy_WT,t,'WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);
    if isempty(W_WT)
        initPKWdth_WT(i)=NaN;
    else
        initPKWdth_WT(i)=W_WT(1);
    end
    
    
end
 
figure(2)
% plot(10.^(logcAMP),initPKWdth./initPKWdth(1),'.','color',[0.1 0.5 0.25],'MarkerSize',15);
plot(10.^(logcAMP),initPKWdth_WT,'.','color',[0.1 0.5 0.25],'MarkerSize',15);
hold on

plot(10.^(logcAMP),initPKWdth_PKAR,'.b','MarkerSize',15);
xlabel('[cAMP]_{ex}')
ylabel('Time (Abs. U.)')
title('PInitial peak width')
set(gca,'Xscale','log','FontSize',10)
% set(gca,'Xscale','log','yLim',[1 1.4],'FontSize',10)
legend('WT','PKAR mutant')

%% %% fig. 2 E Mean Oscilation time
logcAMP=2:0.33:4;
n=numel(logcAMP);
MeanOscPk=zeros(1,n);
t=0:dt:1500;
StimInitTime=dt;
for i=1:n
cAMP=10.^(logcAMP(i));
stim=zeros(size(t));stim(StimInitTime/dt:end)=cAMP;
[tout,xout]=euler_solver(@Sgro_single_cell_new,t,[A0;R0],param_WT,stim);
A_WT=xout(1,:);
yy_PKAR = smooth(A_WT,200,'moving'); % moving average filter
 plot(t,A_WT,'color',[0 1 0],'LineWidth',3); hold on
 %findpeaks(yy,t,'Annotate','extents','WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);hold on
[PKS,LOCS,W]=findpeaks(yy_PKAR,t,'WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);
Osctime=diff(LOCS);
MeanOscPk(i)=mean(Osctime);
end

figure(3)
MeanOscPk=MeanOscPk./MeanOscPk(1);
plot(10.^(logcAMP),MeanOscPk,'.','color',[0.1 0.5 0.25],'MarkerSize',25);
xlabel('[cAMP]_{ex}')
ylabel('Time (Rel. U.)')
set(gca,'Xscale','log','FontSize',20,'Xlim',[0.9e2 1.1e4],'Ylim',[0.38 1.1])
%% Fig 3 not the same with paper
dt=0.005;t=0:dt:30*param_WT.Nt;
cAMP=10000;
stim=zeros(size(t));
param_WT.e=0.1;

tStepOn=floor(1.25./dt.*param_WT.Nt); % at this time step step trigger happens
tStepOff=floor(3.75./dt.*param_WT.Nt);% at this time step step trigger stops
stim(tStepOn:tStepOff)=cAMP;
% % linear ramp
% tRampOn=10./dt.*param_WT.Nt; % at this time step ramp trigger happens
% tRampOff=40./dt.*param_WT.Nt;% at this time step ramp trigger stops
% k=(cAMP-0)./(tRampOff-tRampOn); b=0-k.*tRampOn;
% tRamp=linspace(tRampOn,tRampOff,tRampOff-tRampOn+1); stim(tRampOn:tRampOff)=k*tRamp+b;

% exponential ramp, take the left tail of y=a*exp(bx) and fit function
tRampOn=10./dt.*param_WT.Nt; % at this time step ramp trigger happens
tRampOff=25./dt.*param_WT.Nt;% at this time step ramp trigger stops
tRamp=linspace(tRampOn,tRampOff,tRampOff-tRampOn+1); 

a=log(cAMP+1)/15; % Ramp stimulus at t=0: cAMP=exp(a*t)-1
temp_t=linspace(0,15,length(tRamp));
stim(tRampOn:tRampOff)=exp(a.*temp_t)-1;

tStepOn=25./dt.*param_WT.Nt; % at this time step trigger happens
tStepOff=30./dt.*param_WT.Nt;% at this time step trigger stops
stim(tStepOn:tStepOff)=cAMP;

A0=-1.5;R0=-0.5;
[tout,xout]=euler_solver(@Sgro_single_cell_new,t,[A0;R0],param_PKAR,stim);

% A_WT_orig=xout(1,:);R_WT_orig=xout(2,:); % un-scaled, un-offset A and R after simulation
% A_WT=(A_WT_orig+param_WT.offset_A)./param_WT.Na;
% R_WT=R_WT_orig./param_WT.Na; 

A_PKAR_orig=xout(1,:);R_PKAR_orig=xout(2,:); % un-scaled, un-offset A and R after simulation
A_PKAR=(A_PKAR_orig+param_WT.offset_A)./param_WT.Na;
R_PKAR=R_PKAR_orig./param_WT.Na; 
t_plot=tout./param_WT.Nt;

figure();  
subplot(2,1,1)
plot(t_plot,stim)
%ylim([0 310])
subplot(2,1,2)
plot(t_plot,A_PKAR,'color',[0.1 0.5 0.25],'LineWidth',3); hold on
% line([StimInitLine StimInitLine],get(h,'YLim'),'Color',[0 1 0],'LineWidth',3,'LineStyle','--')
ylabel('Amplitude'); xlabel('Time (T)')
title('step and ramp response to 10000nM cAMP, PKAR')
