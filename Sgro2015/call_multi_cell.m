% Sgro 2015 model. call functions and solvers to reproduce Fig 5- 7
% clear; clc; close all;

% parameters
para.e=0.1; %tauA/tauR: episilon, ratio between activator and repressor time scales
para.tauA= 0.09;
para.tauR=para.tauA/para.e;
para.g=0.5; %gamma, repressor degradation rate
para.c0=1.2; % steady state repressor value in the absense of external cAMP
para.sigma=0.15; %0.15; %sigma, noise strength
para.N=100; % number of cells in population simulations
para.a=0.058; % cAMP response magnitude
para.alpha0=800;% basal cAMP leakage rate
para.alpha_pde= 1000; % basal PDE leakage rate
para.Kd= 1e-5; % cAMP response threshold
para.S=1e6; % firing cell cAMP release rate  
% simulation time vector
dt=0.005;
t = 0:dt:500;  

% Fig 5- population firing rate relates to cAMP out flow rate j and cell density rho

j= 0:0.1:1;%j= 0:0.05:1;%0.205128;%0.05;%% ;% 0.615384615384615;
logrho=-5.5:0.5:-3; %logrho=-5.5:0.1:-3;%-3.07317073170732;% % %% -4.17073170731707;% -5.02439024390244; % %log of 10
alphafval=2.5; % cAMP feed-in rate
% [x,y,z]=meshgrid(j,rho,alphaf);

%alphaf=alphafval.*ones(n,m);
time_separation=1; 
tic
for p=1:length(j) % 3 %5 %
    for q= 1:length(logrho) %21 %11 %
        [A,R,cAMP]=multi_cell (alphafval,j(p),10^(logrho(q)),para,dt,t,time_separation);
%         figure()
%         plot(t(cutoff_time:end),cAMP(cutoff_time:end))
        % [A,R,cAMP]=multi_cell (alphafval,j(p),10^(logrho(q)),para,dt,t);
% rate(q,p)=multi_cell (alphafval,j(p),10^(logrho(q)),para,dt,t);
% fprintf('p= %d of %d and q = %d of %d\n', p,n,q,m);
% J(q,p)=j(p);
% Rho(q,p)=logrho(q);
        % only considering after t=50 (after the first two spikes)
        threshold_pop=-1;
        cutoff_time=100/dt;
        Nt=27; % normalization factor of t
%         threshold=0;
%         for i=1:para.N
%             [pks,locs]=findpeaks(A(i,cutoff_time:end),t(cutoff_time:end),'MinPeakProminence',0.5);
%             NumOfPeaks(i)=sum(pks>threshold);     
%             if mod(i,10)==0
%                 plot(t,A(i,:));hold on;
%             end
%         end
%         NumberOfPeaksMeanSingle=mean(NumOfPeaks);
%         single_rate_mean(p,q)= NumberOfPeaksMeanSingle./(t(end)-cutoff_time*dt).*Nt;    
        A_pop_avg=mean(A); % A_pop_avg=mean(A(:,10001:end));  
        [pks,locs]=findpeaks(A_pop_avg(cutoff_time:end),t(cutoff_time:end),'MinPeakProminence',0.5);    
%         figure()
%         findpeaks(A_pop_avg(cutoff_time:end),t(cutoff_time:end),'MinPeakProminence',0.5);         
        PopNumOfPeaks=sum(pks>threshold_pop);
        pop_rate(q,p)=PopNumOfPeaks./(t(end)-cutoff_time*dt).*Nt;
    end
    disp(strcat('Simulation ended for j=',num2str(j(p))));
    toc
end
%  fig 5A
figure(2)
j_surf=repmat(j,length(logrho),1);
logrho_surf=repmat(logrho',1,length(j));
surf(j_surf,logrho_surf,pop_rate)
% surf(j_surf,logrho_surf,pop_rate(1:9,1:11));
title(['noise=',num2str(para.sigma),',time separation ',num2str(time_separation)])
view(2); colormap(jet); colorbar;caxis([0 0.6]);

%% Save variables to curret folder
filename = 'pop_firing_rate_vars.mat';
save(filename)
%% Fig 5C
rho=10.^(logrho);
figure()
for p=7:length(j)
    hold all
    log_rho_j(p,:)=log10(rho./j(p));
    plot(log_rho_j(p,:),pop_rate(:,p),'.','MarkerSize',15)
%     legend(strcat('J= ',num2str(j(p))))
%     hold all
% holds old plot for multi-line plot
	legendText{p-6} = sprintf(strcat('J= ', num2str(j(p))));
	legend(legendText);
	drawnow; % Force screen refresh.
end
xlim([-6 0]); ylim([0 0.7]);
xlabel('log_{10}(\rho/J)'); ylabel('Population firing rate')

