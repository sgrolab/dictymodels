% Gregor 2010- phase equation model
% Fig 4A
% clear all
% clc
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
% clear all
% clc

t=0:0.01:100;
Nc=50; % number of cells in the population
Amax=20;Abas=0.4;

% log10rho_space= linspace(-4.5,-1,4); rho_space=10.^log10rho_space;
% log10k_space=linspace(0.02,1,4); % ml/min
% k_space=10.^log10k_space;
% log10rho_space= linspace(-3.5,1,26); rho_space=10.^log10rho_space;
log10rho_space= linspace(-3.5,1,26); rho_space=10.^log10rho_space;
% log10k_space=linspace(0.02,2,6); % ml/min
% k_space=10.^log10k_space;
% k_space=linspace(1,100,21);
k_space=linspace(1,100,21);

pop_rate_Gregor=nan(length(log10rho_space),length(k_space));
eta=0.06; % noise strength£¬bigger than stated in the paper
camp_cyt0=Abas*ones(Nc,1);

tic
for p=1:length(k_space)
    for q=1:length(log10rho_space)
        [camp_cyt,camp_ex,thetai,noise]=Gregor2010_fun(camp_cyt0,rho_space(q),k_space(p),t,eta);
        camp_cyt_pop_mean=mean(camp_cyt);
        % check traces
%         camp_cyt_smooth = smooth(camp_cyt_pop_mean,100,'moving'); % moving average filter
%          plot(t,camp_cyt_smooth,'r');hold on
%          figure();plot(t,camp_cyt_pop_mean,'b')
        [PKS,LOCS,W]=findpeaks(camp_cyt_pop_mean,t,'WidthReference','halfheight','MinPeakHeight',1,'MinPeakDistance',1,'MinPeakProminence',0.5);
        PopNumOfPeaks_Gregor=length(PKS);
        pop_rate_Gregor(q,p)=PopNumOfPeaks_Gregor./max(t);
%         if isempty(LOCS) % when there is no oscillation
%             pop_rate_Gregor(q,p)=0;
%         else
%             Osctime=diff(LOCS);Osctime(1:floor(length(Osctime)./2))=[];
%             pop_rate_Gregor(q,p)=mean(Osctime); 
%         end
    end
    disp(strcat('Simulation ended for k=',num2str(k_space(p))));
    toc
end
% Phase plot
figure()
k_surf=repmat(k_space,length(log10rho_space),1);
logrho_Gregor_surf=repmat(log10rho_space',1,length(k_space));
surf(k_surf,logrho_Gregor_surf,pop_rate_Gregor);
title(['noise=',num2str(eta)])
view(2); colormap(jet); colorbar;%caxis([0 0.6]);
    



