% Sgro 2015 single cell
clear all;clc;
% parameters
param.e=0.2; %tauA/tauR; %episilon
param.tauA=0.09;
param.tauR=param.tauA/param.e;
param.g=0.5; %gamma
param.c0=1.2;
param.sigma=0.15; %sigma
param.N=100; % number of cells in population simulations
param.a=0.058;
param.alpha0=800;
param.alpha_pde= 1000;
param.Kd= 1e-5;
param.S=1e6;                          

param.Nt=27; % normalization factor of t
param.Na=3.5;  % normalization factor of A
param.offset_A=1.5;

%% fig 1 G& H
dt=0.005;t=0:dt:6.5*param.Nt;
cAMP=1;
StimInitTime=1; % stimulation time after Nt normalization 
stim=zeros(size(t));stim(StimInitTime/dt*param.Nt:end)=cAMP;
A0=-1.5;R0=-0.5;
[tout,xout]=euler_solver(@Sgro_single_cell_new,t,[A0;R0],param,stim);

A=xout(1,:);R=xout(2,:);
A_orig=A;R_orig=R; % un-scaled, un-offset A and R after simulation
A=(A+param.offset_A)./param.Na;
R=R./param.Na;
t_plot=tout./param.Nt;
StimInitLine=StimInitTime;

figure();  h=axes;
plot(t_plot,A,'color',[0.1 0.5 0.25],'LineWidth',3); hold on
line([StimInitLine StimInitLine],get(h,'YLim'),'Color',[0 1 0],'LineWidth',3,'LineStyle','--')
ylabel('Amplitude'); xlabel('Time (T)')
set(gca,'FontSize',20,'xLim',[0 6.5])

%% Fig 2C- initial peak width
logcAMP=0:0.33:4;
n=numel(logcAMP);
initPKWdth=zeros(1,n);
t=0:dt:200;
A0=-1.5;R0=-0.5;
StimInitTime=1;
for i=1:n% n=1 to test findpeaks 
cAMP=10.^(logcAMP(i));
stim=zeros(size(t));stim(StimInitTime/dt*param.Nt:end)=cAMP;
[tout,xout]=euler_solver(@Sgro_single_cell_new,t,[A0;R0],param,stim);
A=xout(1,:);
yy = smooth(A,200,'moving'); % moving average filter
% plot(t,A,'color',[0 1 0],'LineWidth',3); hold on % test findpeaks
% findpeaks(yy,t,'Annotate','extents','WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);hold on % test findpeaks
[PKS,LOCS,W]=findpeaks(yy,t,'WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);
initPKWdth(i)=W(1);
end
figure(2)
plot(10.^(logcAMP),initPKWdth./initPKWdth(1),'.','color',[0.1 0.5 0.25],'MarkerSize',15);
xlabel('[cAMP]_{ex}')
ylabel('Time (Rel. U.)')
set(gca,'Xscale','log','yLim',[1 1.4],'FontSize',10)

%% %% fig. 2 E Mean Oscilation time
logcAMP=2:0.33:4;
n=numel(logcAMP);
MeanOscPk=zeros(1,n);
t=0:dt:1500;
StimInitTime=dt;
for i=1:n
cAMP=10.^(logcAMP(i));
stim=zeros(size(t));stim(StimInitTime/dt:end)=cAMP;
[tout,xout]=euler_solver(@Sgro_single_cell_new,t,[A0;R0],param,stim);
A=xout(1,:);
yy = smooth(A,200,'moving'); % moving average filter
 plot(t,A,'color',[0 1 0],'LineWidth',3); hold on
 %findpeaks(yy,t,'Annotate','extents','WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);hold on
[PKS,LOCS,W]=findpeaks(yy,t,'WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);
Osctime=diff(LOCS);
MeanOscPk(i)=mean(Osctime);
end

figure(3)
MeanOscPk=MeanOscPk./MeanOscPk(1);
plot(10.^(logcAMP),MeanOscPk,'.','color',[0.1 0.5 0.25],'MarkerSize',25);
xlabel('[cAMP]_{ex}')
ylabel('Time (Rel. U.)')
set(gca,'Xscale','log','FontSize',20,'Xlim',[0.9e2 1.1e4],'Ylim',[0.38 1.1])
%% Fig 3 not really working...... Not the same with paper
dt=0.005;t=0:dt:30*param.Nt;
cAMP=300;
stim=zeros(size(t));
param.e=0.1;

tStepOn=floor(1.25./dt.*param.Nt); % at this time step step trigger happens
tStepOff=floor(3.75./dt.*param.Nt);% at this time step step trigger stops
stim(tStepOn:tStepOff)=cAMP;
% % linear ramp
% tRampOn=10./dt.*param.Nt; % at this time step ramp trigger happens
% tRampOff=25./dt.*param.Nt;% at this time step ramp trigger stops
% k=(cAMP-0)./(tRampOff-tRampOn); b=0-k.*tRampOn;
% tRamp=linspace(tRampOn,tRampOff,tRampOff-tRampOn+1); stim(tRampOn:tRampOff)=k*tRamp+b;

% exponential ramp, take the left tail of y=a*exp(bx) and fit function
tRampOn=10./dt.*param.Nt; % at this time step ramp trigger happens
tRampOff=25./dt.*param.Nt;% at this time step ramp trigger stops
b=log(cAMP./0.001)./(tRampOff-tRampOn); a=cAMP./exp(b.*tRampOff);
tRamp=linspace(tRampOn,tRampOff,tRampOff-tRampOn+1); stim(tRampOn:tRampOff)=a*exp(b.*tRamp);

tStepOn=25./dt.*param.Nt; % at this time step trigger happens
tStepOff=27.5./dt.*param.Nt;% at this time step trigger stops
stim(tStepOn:tStepOff)=cAMP;

A0=-1.5;R0=-0.5;
[tout,xout]=euler_solver(@Sgro_single_cell_new,t,[A0;R0],param,stim);

A=xout(1,:);R=xout(2,:);
A_orig=A;R_orig=R; % un-scaled, un-offset A and R after simulation
A=(A+param.offset_A)./param.Na;
R=R./param.Na;
t_plot=tout./param.Nt;

figure();  
subplot(2,1,1)
plot(t_plot,stim)
ylim([0 310])
subplot(2,1,2)
plot(t_plot,A,'color',[0.1 0.5 0.25],'LineWidth',3); hold on
% line([StimInitLine StimInitLine],get(h,'YLim'),'Color',[0 1 0],'LineWidth',3,'LineStyle','--')
ylabel('Amplitude'); xlabel('Time (T)')
