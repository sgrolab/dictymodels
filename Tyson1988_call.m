%% Tyson 1988
clear all;clc;
% Some parameters based on Table 2, parameter set B
global L1 L2 k c alpha lambda1 lambda2 s1 s2 s epsilon epsilon_prime D k1 ke
% Parameter set B
L1=10; L2=0.005; k=18.5; c=10; alpha=3; lambda1=1e-3; lambda2=2.4; 
s1=950; s2=0.05; s=0.07; epsilon=0.019; epsilon_prime=0.01; D=0.024;% mm2/min
k1=0.036;ke=3.6; 
% set up mesh grid
global X Y x y numofgrids h
[X,Y]=meshgrid(0:0.1:9,0:0.1:9);
x=k1./sqrt(ke*D)*X;
y=k1./sqrt(ke*D)*Y; % dimensionless spatial variables
numofgrids=length(X(1,:));
h=x(1,2)-x(1,1);
%% Set initial conditions on spatial mesh, collapse into one vector
rho0=reshape(0.7*ones(numofgrids),[],1);
beta0=reshape(zeros(numofgrids),[],1);
gamma0_mat=zeros(numofgrids);
% set the area and amplitude of initial [cAMP]external, area is in the
% middle of the mesh grids
a=10;
gamma0_mat((floor(numofgrids/2)-a):(floor(numofgrids/2)+a),(floor(numofgrids/2)-a):(floor(numofgrids/2)+a))=0.2.*ones((2*a+1),(2*a+1)); 
gamma0=reshape(gamma0_mat,[],1);

x0=[rho0;beta0;gamma0];
% set tspan and calculation step size
tspan=0:0.1:80;
% solve using ode45
[t1,x1]=ode45(@(t,x) Tyson1988_fun(t,x),tspan,x0);
%% reconstruct matrix from calculated vectors
[t_length,x1_length]=size(x1);
for i=1:1:t_length
    rho(:,:,i)=reshape(x1(i,1:numofgrids^2),numofgrids,numofgrids);
    beta(:,:,i)=reshape(x1(i,(numofgrids^2+1):numofgrids.^2*2),numofgrids,numofgrids);
    gamma(:,:,i)=reshape(x1(i,(numofgrids.^2*2+1):numofgrids.^2*3),numofgrids,numofgrids);
end
% save variables in workspace
save test052918.m

% Test: surf of beta at one time frame
surf(X,Y,beta(:,:,700));
xlabel('um')
ylabel('um')
zlabel('beta: [cAMP]in/K_R')
title(['Time frame #: ', num2str(i)])

%% Make video: beta change spatial temporally
v = VideoWriter('beta_cAMP_in.avi');
open(v);
surf(X,Y,beta(:,:,1));

set(gca,'nextplot','replacechildren'); 
for i = 400:1:t_length
   surf(X,Y,beta(:,:,i))
   % view(2)
   xlim([4 6])
   ylim([4 6])
   xlabel('um')
   ylabel('um')
   zlabel('beta: [cAMP]in/K_R')
   title(['Beta time frame #: ', num2str(i)])
   colorbar
   frame = getframe;
   writeVideo(v,frame);
end

