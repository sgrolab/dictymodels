clc; clear; close all;

%% fig. 1G&H
%cAmp=1 ;         %Small Stimulus ; cAmp=1 nM  Fig 1G
cAmp=10000;     %Small Stimulus ; cAmp=10 uM  Fig 1H
t_tot=175 ;      %Total simulation time in absolute units
[A,R,t,lineVal]=single_cell(cAmp,t_tot);% 

figure(); 
h=axes;
plot(t,A,'color',[0.1 0.5 0.25],'LineWidth',3); hold on
line([lineVal lineVal],get(h,'YLim'),'Color',[0 1 0],'LineWidth',3,'LineStyle','--')
ylabel('Amplitude')
xlabel('Time (T)')
set(gca,'FontSize',14,'xLim',[0 6.5])

figure();
plot(A,R)
%% fig. 2 C - initial peak width

logcAmp=0:0.33:4;
n=numel(logcAmp);
initPKWdth=zeros(1,n);
t_tot=1000;


for i=1:1 % n=1 to test findpeaks 
cAmp=10.^(logcAmp(i));
%cAmp=10^2; % test findpeaks
[A,R,t,lineVal]=single_cell(cAmp,t_tot);
yy = smooth(A,200,'moving'); % moving average filter
% plot(t,A,'color',[0 1 0],'LineWidth',3); hold on % test findpeaks
% findpeaks(yy,t,'Annotate','extents','WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);hold on % test findpeaks
[PKS,LOCS,W]=findpeaks(yy,t,'WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);
initPKWdth(i)=W(1);
end
figure(2)
plot(10.^(logcAmp),initPKWdth./initPKWdth(1),'.','color',[0.1 0.5 0.25],'MarkerSize',25);
xlabel('[cAMP]_{ex}')
ylabel('Time (Rel. U.)')
set(gca,'Xscale','log','yLim',[1 1.4],'FontSize',14)

%% %% fig. 2 E Mean Oscilation time
logcAmp=2:0.33:4;
n=numel(logcAmp);
MeanOscPk=zeros(1,n);
t_tot=15000;
for i=1:n
cAmp=10.^(logcAmp(i));
%cAmp=10^2;
[A,R,t,lineVal]=single_cell(cAmp,t_tot);
yy = smooth(A,200,'moving'); % moving average filter
% plot(t,A,'color',[0 1 0],'LineWidth',3); hold on
% findpeaks(yy,t,'Annotate','extents','WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);hold on
[PKS,LOCS,W]=findpeaks(yy,t,'WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);
Osctime=diff(LOCS);
MeanOscPk(i)=mean(Osctime);
end

figure(3)
MeanOscPk=MeanOscPk./MeanOscPk(1);
plot(10.^(logcAmp),MeanOscPk,'.','color',[0.1 0.5 0.25],'MarkerSize',25);
xlabel('[cAMP]_{ex}')
ylabel('Time (Rel. U.)')
set(gca,'Xscale','log','FontSize',14,'Xlim',[0.9e2 1.1e4],'Ylim',[0.38 1.1])


