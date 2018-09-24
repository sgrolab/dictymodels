% Gregor 2010- phase equation model
% Fig 4A
clear all
clc
dt=0.1;
k=5;% ml/min
Nc=100; % number of cells in the population
noise=0.002; % noise larger than in the paper

t1=0:dt:50;
camp_cyt0=0.4*ones(Nc,1);
rho1=0;% ML, monolayer
[camp_cyt1,camp_ex1,thetai1,noise1]=Gregor2010_fun(camp_cyt0,rho1,k,t1,noise);

t2=50+dt:dt:83;
rho2=1/768;% ML, monolayer
[camp_cyt2,camp_ex2,thetai2,noise2]=Gregor2010_fun(camp_cyt1(:,end),rho2,k,t2,noise);

t3=83+dt:dt:125;
rho3=1/192;% ML, monolayer
[camp_cyt3,camp_ex3,thetai3,noise3]=Gregor2010_fun(camp_cyt2(:,end),rho3,k,t3,noise);

t4=125+dt:dt:163;
rho4=1/48;% ML, monolayer, different with paper
[camp_cyt4,camp_ex4,thetai4,noise4]=Gregor2010_fun(camp_cyt3(:,end),rho4,k,t4,noise);

t5=163+dt:dt:210;
rho5=1/12;% ML, monolayer
[camp_cyt5,camp_ex5,thetai5,noise5]=Gregor2010_fun(camp_cyt4(:,end),rho5,k,t5,noise);

t6=210+dt:dt:250;
rho6=1/3;% ML, monolayer
[camp_cyt6,camp_ex6,thetai6,noise6]=Gregor2010_fun(camp_cyt5(:,end),rho6,k,t6,noise);

Amax=20;
t=[ t2 t3 t4 t5 t6];
camp_cyt=[camp_cyt2,camp_cyt3,camp_cyt4,camp_cyt5,camp_cyt6 ];
camp_ex=[camp_ex2,camp_ex3,camp_ex4,camp_ex5,camp_ex6];

figure()
[hAx,hLine1,hLine2]=plotyy(t,camp_cyt(1,:),t,camp_ex);hold on
set(hAx(1),'Ylim',[0 21]);
set(hAx(2),'Yscale','log','box','off','Ylim',[0 10]);
set(hLine2,'LineWidth',1.5,'Color',[0 0 0])
title('Fig 4a')
xlabel('Time (min)')
ylabel(hAx(1),'[cAMP]_{cyt}') % left y-axis 
ylabel(hAx(2),'[cAMP]_{ext}') % right y-axis

for i=1:Nc
    plot(t,camp_cyt(i,:),'LineWidth',0.3)
    hold on
end
% dash line showing density change time
line([83 83],[0 22],'Color',[1 0 0],'LineWidth',1,'LineStyle','--'); hold on
line([125 125],[0 22],'Color',[1 0 0],'LineWidth',1,'LineStyle','--'); hold on
line([163 163],[0 22],'Color',[1 0 0],'LineWidth',1,'LineStyle','--'); hold on
line([210 210],[0 22],'Color',[1 0 0],'LineWidth',1,'LineStyle','--'); hold on

%% Fig 4B
clear all
clc

t=0:0.1:100;
Nc=50; % number of cells in the population
Amax=20;Abas=0.4;

log10rho_k= linspace(-3.5,-0.1,30);
rho_k_space=10.^log10rho_k;
k=5; % ml/min
eta=0.02; % noise strength£¬bigger than stated in the paper
rho_space=rho_k_space.*k;

camp_cyt0=Abas*ones(Nc,1);

for i=1:length(rho_k_space)
    % with noise
    [camp_cyt,camp_ex,thetai,noise]=Gregor2010_fun(camp_cyt0,rho_space(i),k,t,eta);
    for j=1:Nc
%         camp_cyt_smooth(j,:) = smooth(camp_cyt(j,:),100,'moving'); % moving average filter
%          plot(t,camp_cyt_smooth(j,:),'r');hold on
%          plot(t,camp_cyt(j,:),'b')
     [PKS,LOCS,W]=findpeaks(camp_cyt(j,:),t,'WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);
   
    if isempty(LOCS) % when there is no oscillation
        MeanT_SC(j,i)=0;
    else
        Osctime=diff(LOCS);Osctime(1:floor(length(Osctime)./2))=[];
        MeanT_SC(j,i)=mean(Osctime); % single cell oscillation period
    end
     end
    MeanT_MC_noise=mean(MeanT_SC); % averaged oscillation period for multiple cells
    
%     % without noise 
    [camp_cyt,camp_ex,thetai,noise]=Gregor2010_fun(camp_cyt0,rho_space(i),k,t,0);% without noise
    for j=1:Nc
%         camp_cyt_smooth(j,:) = smooth(camp_cyt(j,:),100,'moving'); % moving average filter
%          plot(t,camp_cyt_smooth(j,:),'r');hold on
%          plot(t,camp_cyt(j,:),'b')
     [PKS,LOCS,W]=findpeaks(camp_cyt(j,:),t,'WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);
    if isempty(LOCS) % when there is no oscillation
        MeanT_MC(j,i)=0;
    else
        Osctime=diff(LOCS);Osctime(1:floor(length(Osctime)./2))=[];
        MeanT_MC(j,i)=mean(Osctime); % single cell oscillation period
    end
    end
    MeanT_MC_wonoise=mean(MeanT_MC); % averaged oscillation period for multiple cells

end
% convert period to pulse/ min
for i=1:length(MeanT_MC_noise)
    pulse_noise(i)=1./MeanT_MC_noise(i);
    pulse_noise(isnan(pulse_noise))=0;
    pulse_wonoise(i)=1./MeanT_MC_wonoise(i);
    pulse_wonoise(isinf(pulse_wonoise))=0;
end
% plot figure 4b
plot(log10rho_k, pulse_noise,'b.','MarkerSize',10);hold on
plot(log10rho_k, pulse_wonoise,'kd','MarkerSize',6);
legend(strcat('with noise \eta=',num2str(eta)),'without noise','Location','northwest')
xlabel('log_{10} \rho/k')
ylabel('cAMP pulse/min')
title('Fig 4B')


