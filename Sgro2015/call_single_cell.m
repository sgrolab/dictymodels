clc; clear; close all;
% fig. 1G&H
% cAMP=1 ;         %Small Stimulus ; cAmp=1 nM  Fig 1G
cAMP=10000;     %Small Stimulus ; cAmp=10 uM  Fig 1H
t_tot=175 ;      %Total simulation time in absolute units
[A,R,t,lineVal,A_orig,R_orig]=single_cell_new(cAMP,t_tot);% 


figure(); 
h=axes;
plot(t,A,'color',[0.1 0.5 0.25],'LineWidth',3); hold on
line([lineVal lineVal],get(h,'YLim'),'Color',[0 1 0],'LineWidth',3,'LineStyle','--')
ylabel('Amplitude')
xlabel('Time (T)')
set(gca,'FontSize',20,'xLim',[0 6.5])

% nullcline 
A_null=linspace(-2.5,2.5,200);
dAdt_null_no_stim=A_null-1/3*A_null.^3;% dA/dt=0 nullcline w/o stim
a=0.058;Kd= 1e-5;
dAdt_null_stim=A_null-1/3*A_null.^3+a*log(1+cAMP/Kd)*ones(size(A_null));% dA/dt=0 nullcline with stim
c0=1.2;g=0.5; %gamma
dRdt_null=1/g*(A_null+c0*ones(size(A_null))); % dR/dt=0 nullcline

figure();
plot(A_null,dRdt_null,'r','LineWidth',2)
hold on
plot(A_null,dAdt_null_no_stim,'g','LineWidth',2)
hold on 
plot(A_null,dAdt_null_stim,'k','LineWidth',2)
hold on
plot(A_orig,R_orig,'k','LineWidth',2)
ylabel('Repressor (R)')
xlabel('Activator (A)')
set(gca,'yLim',[-1 2.5],'FontSize',10)
%%
A_null=linspace(-2.5,2.5,200);
dAdt_null_no_stim=A_null-1/3*A_null.^3;% dA/dt=0 nullcline w/o stim
a=0.058;Kd= 1e-5;
dAdt_null_stim=A_null-1/3*A_null.^3+a*log(1+cAMP/Kd)*ones(size(A_null));% dA/dt=0 nullcline with stim
c0=1.2;g=0.5; %gamma
dRdt_null=1/g*(A_null+c0*ones(size(A_null))); % dR/dt=0 nullcline

figure()
subplot(1,2,1)
plot(A_null,dRdt_null,'r','LineWidth',2)
hold on
plot(A_null,dAdt_null_no_stim,'g','LineWidth',2)
hold on 
plot(A_null,dAdt_null_stim,'k','LineWidth',2)
hold on

for i=1:length(t)
    if mod(i,200)==0
        subplot(1,2,1)
        scatter(A_orig(i),R_orig(i),'b.')
        hold on
        
        subplot(1,2,2)
        % h=axes;
        scatter(t(i),A(i)); hold on
        % line([lineVal lineVal],get(h,'YLim'),'Color',[0 1 0],'LineWidth',3,'LineStyle','--');hold on
        ylabel('Amplitude')
        xlabel('Time (T)')
        set(gca,'FontSize',20,'xLim',[0 6.5])
        
        pause(0.0001)
    end
end
title(['\epsilon=',])
%% fig. 2 C - initial peak width

logcAMP=0:0.33:4;
n=numel(logcAMP);
initPKWdth=zeros(1,n);
t_tot=1000;


for i=1:n% n=1 to test findpeaks 
cAMP=10.^(logcAMP(i));
% cAmp=10^2; % test findpeaks
[A,R,t,lineVal,A_orig,R_orig]=single_cell_new(cAMP,t_tot);
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
t_tot=15000;
for i=1:n
cAMP=10.^(logcAMP(i));
%cAmp=10^2;
[A,R,t,lineVal]=single_cell_old(cAMP,t_tot);
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

%% %% Mean Oscilation height
logcAMP=2:0.33:4;
n=numel(logcAMP);
MeanOscPk=zeros(1,n);
t_tot=1000;
for i=1:n
    cAMP=10.^(logcAMP(i));
    %cAmp=10^2;
    [A,R,t,lineVal,A_orig,R_orig]=single_cell_old(cAMP,t_tot);
    A=A_orig;
    yy = smooth(A,200,'moving'); % moving average filter
    % plot(t,A,'color',[0 1 0],'LineWidth',3); hold on
    % findpeaks(yy,t,'Annotate','extents','WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);hold on
    [PKS,LOCS,W,PROMS]=findpeaks(yy,t,'WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);

    MeanHeight(i)=mean(PROMS);% mean peak hight
    SEMHeight(i)=std(PROMS)./sqrt(length(PROMS));
end
figure(3)
errorbar(10.^(logcAMP), MeanHeight, SEMHeight,'LineWidth',2)
grid
xlabel('[cAMP]_{ex}')
ylabel('1st cAMP spike height (Rel. U.)')
set(gca,'Xscale','log','FontSize',10,'Xlim',[0.9e2 1.1e4],'Ylim',[3.7 3.82])

%% %% fig 3 
cAMP=1 ;         %Small Stimulus ; cAmp=1 nM  Fig 3A
% cAMP=10000;     %Small Stimulus ; cAmp=10 uM  Fig 3B
t_tot=540 ;      %Total simulation time in absolute units
[A,R,t,lineVal,A_orig,R_orig]=single_cell_new(cAMP,t_tot);% 

figure(); 
h=axes;
plot(t,A,'color',[0.1 0.5 0.25],'LineWidth',3); hold on
% line([lineVal lineVal],get(h,'YLim'),'Color',[0 1 0],'LineWidth',3,'LineStyle','--')
ylabel('Amplitude')
xlabel('Time (T)')
set(gca,'FontSize',20,'xLim',[0 20])
%% nullcline 
A_null=linspace(-2.5,2.5,200);
dAdt_null_no_stim=A_null-1/3*A_null.^3;% dA/dt=0 nullcline w/o stim
a=0.058;Kd= 1e-5;
dAdt_null_stim=A_null-1/3*A_null.^3+a*log(1+cAMP/Kd)*ones(size(A_null));% dA/dt=0 nullcline with stim
c0=1.2;g=0.5; %gamma
dRdt_null=1/g*(A_null+c0*ones(size(A_null))); % dR/dt=0 nullcline

Nt=27; dt=0.005; tTrig=2.5./dt.*Nt; % at this time step trigger happens
% fig. 3 step increase vs ramp increase
tStepOff=tTrig+Nt.*5./dt;% at this time step trigger stops
tRampOn=tStepOff+Nt.*5./dt; % at this time step ramp trigger happens
tRampOff=tRampOn+Nt.*5./dt;
%%
figure();
plot(A_null,dRdt_null,'r','LineWidth',2)
hold on
plot(A_null,dAdt_null_no_stim,'g','LineWidth',2)
hold on 
plot(A_null,dAdt_null_stim,'k','LineWidth',2)
hold on
plot(A_orig(tRampOn:tRampOff),R_orig(tRampOn:tRampOff),'k','LineWidth',2)
ylabel('Repressor (R)')
xlabel('Activator (A)')
set(gca,'yLim',[-1 2.5],'xlim',[-2.5 2.5],'FontSize',10)