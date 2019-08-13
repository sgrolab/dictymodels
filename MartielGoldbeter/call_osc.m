
clc;close all;clear

tspan = [0 34];    
%init values
y1_0 = 0.6;     %rho
y2_0 = 3;       %alpha
y3_0 = 0.9;     %beta
y4_0 = 0;       %gamma

[T,Y] = ode45(@osc,tspan,[y1_0 y2_0 y3_0 y4_0]);    %solve odes

figure(2)
yyaxis right
plot(T,Y(:,1),'r-','linewidth',3) ; hold on
ylim([0 1])
ylabel('\rho_{T}')
yyaxis left
plot(T,Y(:,2),'m-','linewidth',3)


plot(T,Y(:,3)/60,'b-','linewidth',3);hold on
plot(T,Y(:,4)/60,'k-','linewidth',3)
ylim([0 5])
xlim([0 33])
ylabel('\alpha, \beta/60, \gamma/60')
xlabel('Time (min)')
legend({'\alpha','\beta','\gamma','\rho_{T}'})
set(gca,'fontsize',18)
%% %%%%%%%%%%%%plot%%%%%%%%%%%%%%%

subplot (1,2,1); 
yyaxis right
plot(T,Y(:,1),'r','linewidth',3) ; hold on
ylabel('\rho_{T}')
set(gca,'fontsize',18)
xlim([0 33])

subplot (1,2,1); 
ylim([0 1])
yyaxis left
plot(T,Y(:,2),'m','linewidth',3)
ylim([0 5])
xlim([0 33])
legend('\alpha','\rho_{T}')
xlabel('Time (min)')
ylabel('\alpha')
set(gca,'fontsize',18)

subplot (1,2,2);plot(T,Y(:,3),'b','linewidth',3);hold on
subplot (1,2,2);plot(T,Y(:,4),'k','linewidth',3)
xlim([0 34])
legend('\beta','\gamma')
xlabel('Time (min)')
ylabel('\beta, \gamma')
set(gca,'fontsize',18)

figure(2)
yyaxis right
plot(T,Y(:,1),'r-','linewidth',3) ; hold on
ylim([0 1])
ylabel('\rho_{T}')
yyaxis left
plot(T,Y(:,2),'m-','linewidth',3)


plot(T,Y(:,3)/60,'b-','linewidth',3);hold on
plot(T,Y(:,4)/60,'k-','linewidth',3)
ylim([0 5])
xlim([0 33])
ylabel('\alpha, \beta/60, \gamma/60')
xlabel('Time (min)')
legend({'\alpha','\beta','\gamma','\rho_{T}'})
set(gca,'fontsize',18)
%% 
%%%%% Three variable model with parameter input%%%%%
% Some parameters based on Table 2
alpha=3;% intracellular ATP level 
sigma=0.6; 
q=4000; 
ki=1.7; 
kt=0.9;
h=5; 
ke=5.4; 
k1=0.036; 
k2=0.666; 
L1=10;
L2=0.005; 
c=10; 
theta=0.01; % not available from experimental data
lambda=0.01;% not available from experimental data
epsilon=1; % not available from experimental data

x0=[0.7 0 0];
tspan=[0 80];

[t1,x1]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan,x0);

% plot figure 2a
plot(t1,x1(:,1),'k',t1,x1(:,2),'r',t1,x1(:,3),'b')
title('Three variable simulation')
legend('rho_T: total fraction of active receptor','[cAMP]in/K_R','[cAMP]ex/K_R')

% Plot figure 3a
rho_T=x1(:,1);beta=x1(:,2);gamma=x1(:,3);
c=100;k1=0.4; k2=0.04; L1=10;L2=0.1; 

Y_bar=rho_T.*gamma./(ones(size(gamma))+gamma)+(ones(size(gamma))-rho_T).*c.*gamma./((ones(size(gamma))+c.*gamma));
subplot(2,1,1)
plot(t1,Y_bar,'b',t1,rho_T,'r',t1,gamma,'g')
legend('Y_bar: total binding of cAMP to receptor R and D','rho_T: fraction of in active receptor ','gamma: [cAMP]ex/KR')

% plot Net modification fluxes
phi_1=-rho_T.*k1./(ones(size(gamma))+gamma)+(ones(size(gamma))-rho_T).*k1*L1./((ones(size(gamma))+c.*gamma));
phi_2=-rho_T.*k2.*gamma./(ones(size(gamma))+gamma)+(ones(size(gamma))-rho_T).*k2*L2.*c.*gamma./((ones(size(gamma))+c.*gamma));
subplot(2,1,1)
plot(t1,phi_1,t1,phi_2)
legend('Net modification fluxes min-1, phi_1','phi_2')
