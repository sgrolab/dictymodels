
clc;close all;clear
% Sustained oscillations, 4 variable model, fig. 2
tspan = [0 34];    
%init values
y1_0 = 0.6;     %rho
y2_0 = 3;       %alpha
y3_0 = 0.9;     %beta
y4_0 = 0;       %gamma

[T,Y] = ode45(@osc,tspan,[y1_0 y2_0 y3_0 y4_0]);    %solve odes

%%%%%%%%%%%%%plot%%%%%%%%%%%%%%%

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
% Sustained oscillations, 3 variable model, fig. 3
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
tspan=[0 40];

[t1,x1]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan,x0);

% plot figure 2a
figure(3)
plot(t1,x1(:,1),'k',t1,x1(:,2),'r',t1,x1(:,3),'b')
title('Three variable simulation')
legend('rho_T: total fraction of active receptor','[cAMP]in/K_R','[cAMP]ex/K_R')

% figure 3a
figure(4)
rho_T=x1(:,1);beta=x1(:,2);gamma=x1(:,3);

Y_bar=rho_T.*gamma./(ones(size(gamma))+gamma)+(ones(size(gamma))-rho_T).*c.*gamma./((ones(size(gamma))+c.*gamma));
subplot(2,1,1)
plot(t1,Y_bar,'b',t1,rho_T,'r',t1,gamma,'g')
legend('Y_bar: total binding of cAMP to receptor R and D','rho_T: fraction of in active receptor ','gamma: [cAMP]ex/KR')

% plot Net modification fluxes
phi_1=-rho_T.*k1./(ones(size(gamma))+gamma)+(ones(size(gamma))-rho_T).*k1*L1./((ones(size(gamma))+c.*gamma));
phi_2=-rho_T.*k2.*gamma./(ones(size(gamma))+gamma)+(ones(size(gamma))-rho_T).*k2*L2.*c.*gamma./((ones(size(gamma))+c.*gamma));
subplot(2,1,2)
plot(t1,phi_1,t1,phi_2)
legend('Net modification fluxes min-1, phi_1','phi_2')

% fig. 3 b
% Some parameters based on fig. 3b legend
c=100; 
L2=0.1;
k1=0.4;
k2=0.004;
% initial conditions
x0=[0.25 6 0.3];
tspan=[0 40];
% solve differential equations
[t2,x2]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan,x0);
% Plot figure 3b
figure(5)
rho_T=x2(:,1);beta=x2(:,2);gamma=x2(:,3);
Y_bar=rho_T.*gamma./(ones(size(gamma))+gamma)+(ones(size(gamma))-rho_T).*c.*gamma./((ones(size(gamma))+c.*gamma));
subplot(2,1,1)
plot(t2,Y_bar,'b',t2,rho_T,'r',t2,gamma,'g')
legend('Y_bar: total binding of cAMP to receptor R and D','rho_T: fraction of in active receptor ','gamma: [cAMP]ex/KR')

phi_1=-rho_T.*k1./(ones(size(gamma))+gamma)+(ones(size(gamma))-rho_T).*k1*L1./((ones(size(gamma))+c.*gamma));
phi_2=-rho_T.*k2.*gamma./(ones(size(gamma))+gamma)+(ones(size(gamma))-rho_T).*k2*L2.*c.*gamma./((ones(size(gamma))+c.*gamma));
subplot(2,1,2)
plot(t2,phi_1,t2,phi_2)
legend('Net modification fluxes min-1, phi_1','phi_2')

%% Excitability
% Some parameters based on Table 2 and fig 5 legend
alpha=3;% intracellular ATP level 
sigma=0.57; % sigma=0.6; 
q=4000; 
ki=0.958; % ki=1.7; 
kt=0.9;
h=5; 
ke=3.58; % ke=5.4; 
k1=0.036; 
k2=0.666; 
L1=10;
L2=0.005; 
c=10; 
theta=0.01; 
lambda=0.01;
epsilon=0.108; % epsilon=1; 

% fig 5
x0_1=[0.8 0 0];
tspan_1=[-5 0];
[t1,x1]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan_1,x0_1);
x0_2=x1(end,:);x0_2(1,3)=0.3;
tspan_2=[tspan_1(end) 10];
[t2,x2]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan_2,x0_2);
x0_3=x2(end,:);x0_3(1,3)=0.3;
tspan_3=[tspan_2(end) 17.5];
[t3,x3]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan_3,x0_3);
x0_4=x3(end,:);x0_4(1,3)=0.3;
tspan_4=[tspan_3(end) 40];
[t4,x4]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan_4,x0_4);
x_plot=[x1;x2;x3;x4];
t_plot=[t1;t2;t3;t4];

figure(6)
yyaxis right
plot(t_plot,x_plot(:,1),'b','linewidth',1.5) ; hold on
ylabel('\rho_{T}')
set(gca,'fontsize',10)
xlim([-5 40])
ylim([0 1])

yyaxis left
plot(t_plot,x_plot(:,2),'k','linewidth',1.5)
ylim([0 45])
xlabel('Time (min)')
ylabel('\beta')
legend('\rho_{T}','\beta')
title('Fig. 5 Relay response of cAMP pulses')
%% fig 6- threshold for relay, parameters same as fig 5
gamma_plot= linspace(1,22,50);
gamma_space=0.0284.*gamma_plot;% external cAMP stimulation space
beta0=0.03;rho_T0=0.45; % initial conditions can be tuned to match with fig.6
for i=1:length(gamma_space)
    x0_1=[rho_T0 beta0 gamma_space(i)];
    tspan_1=[0 20];
    [t1,x1]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan_1,x0_1);
    Y_bar1=x1(:,1).*x1(:,3)./(ones(size(x1(:,1)))+x1(:,3))+(ones(size(x1(:,1)))-x1(:,1)).*c.*x1(:,3)./((ones(size(x1(:,1)))+c.*x1(:,3)));
    beta=x1(:,2);
%     x0_2=x1(end,:);x0_2(1,3)=0.3;
%     tspan_2=[tspan_1(end) 10];
%     [t2,x2]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan_2,x0_2);
%     Y_bar2=x2(:,1).*x2(:,3)./(ones(size(x2(:,1)))+x2(:,3))+(ones(size(x2(:,1)))-x2(:,1)).*c.*x2(:,3)./((ones(size(x2(:,1)))+c.*x2(:,3)));
    Y_bar_initial(i)=Y_bar1(1);
    [pks,locs,widths,proms] = findpeaks(x1(:,2));
    beta_peak(i)=beta(locs(1))./beta0;
end
% plot(t1,x1(:,2)./beta0)
% (t1,Y_bar1)
figure(7)
yyaxis right
plot(gamma_plot,Y_bar_initial,'b','linewidth',1.5) ; hold on
ylabel('Initial receptor saturation')
set(gca,'fontsize',10)
xlim([1 22])
ylim([0 1.25])

yyaxis left
plot(gamma_plot,beta_peak,'k','linewidth',1.5)
ylim([0 90])
xlabel('magnitude of cAMP pulse \gamma_i/\gamma_0')
ylabel('magnitude of intracellular cAMP response \beta_{max}/\beta_0')
legend('\beta_{max}/\beta_0','Y_bar')
title('Fig. 6 Threshold for relay response')
%% fig 7- ARP and RRP, parameters same as fig.5
delta_t=linspace(0.1,20,50);

x0_1=[0.8 0 0];
tspan_1=[-5 0];
[t1,x1]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan_1,x0_1);
x0_2=x1(end,:);x0_2(1,3)=0.3;
beta1=x1(:,2);
% without the second stim
tspan_2_0=[tspan_1(end) 40];
[t2_0,x2_0]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan_2_0,x0_2);
beta2_0=x2_0(:,2);
[pks1_0,locs1_0,widths1_0,proms1_0] = findpeaks(beta2_0);
loc_1st_peak=locs1_0(1); % Time point when reaches the 1st peak
time_1st_peak=t2_0(loc_1st_peak);
rho_T_0=x2_0(:,1);
% plot(t2_0,rho_T_0)
% hold on 
% plot(t2_0,beta2_0)
beta_1st_peak=beta2_0(locs1_0(1)); % amplitude of the first peak

for i=1:length(delta_t)
    % first peak 
    tspan_2=[tspan_1(end) delta_t(i)];
    [t2,x2]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan_2,x0_2);
    beta2=x2(:,2);
    % second peak
    x0_3=x2(end,:);x0_3(1,3)=0.3;
    tspan_3=[tspan_2(end) delta_t(i)+20];
    [t3,x3]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan_3,x0_3);
    beta3=x3(:,2);
    [pks2,locs2,widths2,proms2] = findpeaks(beta3);
    
   beta_peak_ratio(i)=beta_2nd_peak(i)./beta_1st_peak;
   
   if delta_t(i)<=time_1st_peak
        beta_peak_ratio(i)=0; % for stimulous before the 1st peak is reached, there is absolute refractory period
    end
    
%     else
%         beta_2nd_peak(i)=beta3(locs2(1));
%     end
%     
%     if beta_peak_ratio(i)>1
%         beta_peak_ratio(i)=0;
%     end
end
% plot to test previous code
% plot(t1,beta1,'b')
% hold on
% plot(t2,beta2,'g')
% hold on
% plot(t3,beta3,'r')

figure(7)
yyaxis right
plot(t2_0,rho_T_0,'b','linewidth',1.5) ; hold on
ylabel('\rho_{T}')
set(gca,'fontsize',9)
xlim([0 20])
ylim([0 1.3])

yyaxis left
plot(delta_t,beta_peak_ratio,'k','linewidth',1.5)
ylim([0 1.3])
xlabel('Time interval between simulations, \delta t (min)')
ylabel('Ratio of second to first response \beta_{m2}/ \beta_{m1}')
legend('\beta_{m2}/ \beta_{m1}','\rho_T')
title('Fig. 7 Refractoriness in the relay response of cAMP pulses')
%% Fig 8- 2 variable system vs. 3 variable system
% parameters taken from table 1/ fig. 3
alpha=3;
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
theta=0.01; 
lambda=0.01;
epsilon=1; 
beta_2var=100; % Not specified in the paper
% Initial conditions according to fig 8 caption and try fit with fig 8
x0_3var=[0.6 0.75 0];  
x0_2var=[0.3 0];
tspan=[0 30];

[t1,x1]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan,x0_3var);
[t2,x2]=ode45(@(t,x) two_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon,beta_2var),tspan,x0_2var);

% figure 8
figure
rho_T_3var=x1(:,1);beta_3var=x1(:,2);gamma_3var=x1(:,3);
rho_T_2var=x2(:,1);gamma_2var=x2(:,2);

yyaxis left
plot(t1,gamma_3var,'k--','linewidth',1.5) ; hold on
plot(t2,gamma_2var,'b-','linewidth',1.5) ; hold on
ylabel('\gamma')
set(gca,'fontsize',10)
xlim([0 30])
ylim([0 15])

yyaxis right
plot(t1,rho_T_3var,'g--','linewidth',1.5); hold on
plot(t2,rho_T_2var,'r-','linewidth',1.5);
ylim([0 1])
xlabel('Time (min)')
ylabel('\rho_T')
legend('3 variable \gamma','2 variable \gamma','3 variable \rho_T','2 variable \rho_T')
title('Fig. 6 Comparison of oscillation in 3 and 2 variable systems')

%% Fig 11- Adaptation to constant stimuli
% parameters taken from table 1/ fig. 3
alpha=3;
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
theta=0.01; 
lambda=0.01;
epsilon=1; 
gamma_stim1=0.1; gamma_stim2=1; gamma_stim3=10;
% Initial conditions according to fig 8 caption and try fit with fig 8
% before stimulation
x0_0=[0.8 0 0];
tspan_0=[-1 0];
[t_minus,x_minus]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan_0,x0_0);
% Initial conditions at time 0 when stimulation is applied
x0 =x_minus(end,:);x0(3)=[];
tspan=[0 8];

[t1,x1]=ode45(@(t,x) const_stim(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon,gamma_stim1),tspan,x0);
[t2,x2]=ode45(@(t,x) const_stim(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon,gamma_stim2),tspan,x0);
[t3,x3]=ode45(@(t,x) const_stim(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon,gamma_stim3),tspan,x0);

% figure 11
figure
rho_T_0=x_minus(:,1);beta_0=x_minus(:,2);
rho_T_1=x1(:,1);beta_1=x1(:,2);
rho_T_2=x2(:,1);beta_2=x2(:,2);
rho_T_3=x3(:,1);beta_3=x3(:,2);

subplot(1,2,1)
plot(t_minus,beta_0,'linewidth',1.5) ; hold on
plot(t1,beta_1,'linewidth',1.5) ; hold on
plot(t2,beta_2,'linewidth',1.5) ; hold on
plot(t3,beta_3,'linewidth',1.5) ;
xlabel('Time (min)')
ylabel('Intracellular cAMP, \beta')
xlim([-1 8])
ylim([0 350])
legend('Before cAMP stim', 'a, \gamma=0.1','b, \gamma=1','c, \gamma=10')
title('Fig. 11 Adaptation to constant stimuli')

subplot(1,2,2)
plot(t_minus,ones(size(rho_T_0))-rho_T_0,'linewidth',1.5) ; hold on
plot(t1,ones(size(rho_T_1))-rho_T_1,'linewidth',1.5) ; hold on
plot(t2,ones(size(rho_T_2))-rho_T_2,'linewidth',1.5) ; hold on
plot(t3,ones(size(rho_T_3))-rho_T_3,'linewidth',1.5) ;
xlabel('Time (min)')
ylabel('Intracellular cAMP, \delta_T')
xlim([-1 8])
ylim([0 1])
legend('Before cAMP stim','a, \gamma=0.1','b, \gamma=1','c, \gamma=10')

%% Fig 12. DR curve for constant stim
% parameter values as in fig 11/ fig 3a
gamma_stim_vec=linspace(0,10,50);
% before stimulation
x0_0=[0.8 0 0];
tspan_0=[-1 0];
[t_minus,x_minus]=ode45(@(t,x) three_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon),tspan_0,x0_0);
% Initial conditions at time 0 when stimulation is applied
x0 =x_minus(end,:);x0(3)=[];
tspan=[0 8];
rho_T_0=x_minus(:,1); beta_0=x_minus(:,2);  gamma_0=x_minus(:,3);
Y_bar_0=rho_T_0.*gamma_0./(ones(size(gamma_0))+gamma_0)+ (ones(size(gamma_0))-rho_T_0).*c.*gamma_0./(ones(size(gamma_0))+c.*gamma_0);

for i=1:length(gamma_stim_vec)
    gamma_stim=gamma_stim_vec(i);
    [t1,x1]=ode45(@(t,x) const_stim(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon,gamma_stim),tspan,x0);
    rho_T_1=x1(:,1);
    beta_1=x1(:,2);
    Y_bar_1=rho_T_1.*gamma_stim./(1+gamma_stim)+ (ones(size(beta_1))-rho_T_1).*c.*gamma_stim./(1+c.*gamma_stim);
    Y_bar_1_initial(i)=Y_bar_1(1);% total initial receptor occpation ratio after stimulation
    [pks1,locs1,widths1,proms1] = findpeaks(beta_1);
    if isempty(locs1) % for results without beta peak, set beta maximum to 0
        beta_m(i)=0;
    else
        loc_1st_peak(i)=locs1(1);
        beta_m(i)=beta_1(loc_1st_peak(i)); % height of first intracellular cAMP mpeak
    end
end
% plot to test previous code with one gamma_stim
% plot(t1,beta_1,'b')
% hold on

figure
yyaxis right
plot(gamma_stim_vec,Y_bar_1_initial,'b','linewidth',1.5) ; hold on
ylabel('Y_bar')
set(gca,'fontsize',9)
xlim([0 10])
ylim([0 1.2])

yyaxis left
plot(gamma_stim_vec,beta_m,'k','linewidth',1.5)
ylim([0 400])
xlabel('Constant level of extracellular cAMP, \gamma')
ylabel('Maximum of cAMP response, \beta_{m}')
legend('\beta_{m}','Y_bar')
title('Fig. 12 Dose response curve for constant cAMP stimuli')
%% Fig 13 continuous activation by serial increments in cAMP stim
% parameters the same with fig 11/ 3a
figure
for i=1:25 % 25 steps increase of cAMP stim
    gamma_ex=1.192*10^(-12)*10^7*2^i;
    if i==1
        tspan=[0 1.5];
        x0=[0.8 0];
        beta_c_prev=0.241;
    else
         tspan=[t(end) 1.5*i]; 
         x0=[rho_T(end) beta(end)];
         beta_c_prev=beta_c(end);
    end
    [t,x]=ode45(@(t,x) const_stim(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon,gamma_ex),tspan,x0);
    rho_T=x(:,1); beta=x(:,2); gamma=ones(size(t)).*gamma_ex;
    beta_c=[];beta_c(1,1)=beta_c_prev;
    for j=2:1:length(t)
        beta_c(j,1)=beta_c(j-1,1)+kt*(beta(j)-0.241).*(t(j)-t(j-1));
    end
    % plot the dynamics
    yyaxis left
    plot(t,log(gamma),'k-','linewidth',1.5) ; hold on
    set(gca,'fontsize',10)

    yyaxis right
    plot(t,beta,'r-','linewidth',1.5); hold on
    plot(t,beta_c,'b-','linewidth',1.5);hold on
end
% continuous cAMP stim 
tspan=[t(end) 60];
x0=[rho_T(end) beta(end)];
beta_c_prev=beta_c(end);
[t,x]=ode45(@(t,x) const_stim(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,epsilon,gamma_ex),tspan,x0);
rho_T=x(:,1); beta=x(:,2); gamma=ones(size(t)).*gamma_ex;
beta_c=[];beta_c(1,1)=beta_c_prev;
for j=2:1:length(t)
    beta_c(j,1)=beta_c(j-1,1)+kt*(beta(j)-0.241).*(t(j)-t(j-1));
end
yyaxis left
h1=plot(t,log(gamma),'k-','linewidth',1.5) ; hold on
ylabel('log \gamma')
set(gca,'fontsize',10)
xlim([0 60])
ylim([-11 7])

yyaxis right
h2=plot(t,beta,'r-','linewidth',1.5); hold on
h3=plot(t,beta_c,'b-','linewidth',1.5);hold on
ylim([0 500])
xlabel('Time (min)')
ylabel('Intracellular cAMP, \beta / cAMP secretion, \beta_c')
legend([h1 h2 h3],'\gamma','\beta','\beta_c')
title('Fig. 13 Ccontinuous activation by serial increments in cAMP stim')

