Tyson1988_% Tyson 1988 2 component model.

clear all, clc, close all

add_extra=false; % extra fun

% Some parameters based on Table 2, parameter set B
global L1 L2 k c alpha lambda1 lambda2 s1 s2 s epsilon epsilon_prime D k1 ke
% Parameter set B
L1=10; L2=0.005; k=18.5; c=10; alpha=3; lambda1=1e-3; lambda2=2.4; 
s1=950; s2=0.05; s=0.07; epsilon_prime=0.019; epsilon=0.01; D=0.024;% mm2/min
k1=0.036;ke=3.6; 

% Initializations
global X Y x y numofgrids h
 [X,Y]=meshgrid(0:0.1:9,0:0.1:9); %% field of observation 9mm X 9mm
%  x=k1./sqrt(ke*D)*X;
%  y=k1./sqrt(ke*D)*Y; % dimensionless spatial variables
h=X(1,2)-X(1,1);
numofgrids=length(X(1,:));
% Time step size
dt = 0.01;
% Time vector
t= 0:dt:50;

% Initial position and concentrations for gamma and rho

rho = 0.7.*ones(numofgrids,numofgrids);
rho_new=zeros(numofgrids,numofgrids);
beta = zeros(numofgrids,numofgrids);
beta_new=zeros(numofgrids,numofgrids);
gamma = zeros(numofgrids,numofgrids);
gamma_new = zeros(numofgrids,numofgrids);

% extracellular cAMP (gamma) local boost
w=3;l=3;% witdth and length of the initial [cAMP]ex area
x1=floor(numofgrids/2);x2=floor(numofgrids/2)+w;
y1=floor(numofgrids/2);y2=floor(numofgrids/2)+l;
gamma(y1:y2,x1:x2)=100.*ones((l+1),(w+1)); 

f1 = figure(1)
       subplot(3,1,1)
       imagesc(rho(2:numofgrids-1,2:numofgrids-1))
       title('Initial condition of rho','FontSize', 11, 'FontWeight', 'bold')
       subplot(3,1,2)        
       imagesc(beta(2:numofgrids-1,2:numofgrids-1))
       title('Initial condition of beta','FontSize', 11, 'FontWeight', 'bold')
       subplot(3,1,3)        
       imagesc(gamma(2:numofgrids-1,2:numofgrids-1))
       title('Initial condition of gamma','FontSize', 11, 'FontWeight', 'bold')
movegui(f1,'northwest')

%% Starting simulation

% Storing the f(u,v) and g(u,v)
f1gamma=@(gamma,rho) (1+k.*gamma)./(1+gamma);
f2gamma=@(gamma,rho) (L1+k.*L2.*c.*gamma)./(1+c.*gamma);
phi=@(gamma,rho) (lambda1+(rho.*gamma./(1+gamma)).^2)./(lambda2+(rho.*gamma./(1+gamma)).^2);

for i = 1:length(t)

   % Setting the zero-flux boundaries
   gamma(1,1:numofgrids) = gamma(2,1:numofgrids);     
   gamma(numofgrids,1:numofgrids) = gamma(numofgrids-1,1:numofgrids);
   gamma(1:numofgrids,numofgrids) = gamma(1:numofgrids,numofgrids-1);   
   gamma(1:numofgrids,1) = gamma(1:numofgrids,2);
     
   % Updating u and v with Euler forward and laplacian discretization
   L=4*del2(gamma,h);
   temp=zeros(size(L));
   temp(2:(numofgrids-1),2:(numofgrids-1))=ones((numofgrids-2),(numofgrids-2));
   L=L.*temp;
   
   rho_new= rho+dt.*k1.*(-f1gamma(gamma,rho).*rho+f2gamma(gamma,rho).*(ones(numofgrids,numofgrids)-rho));
   beta_new=beta+dt.*k1.*1./epsilon_prime.*(s1.*phi(gamma,rho)-beta);
   gamma_new = gamma+dt.*(D.*L+((ke.*s.*phi(gamma,rho)-gamma)));
   rho=rho_new;
   beta=beta_new;
   gamma=gamma_new;
   
   
   
   % Plotting every tenth time step
   if (0 == mod(i,5) )       
       hold on
       
       % Plotting gamma and rho
       f2 = figure(2);
       hold on
       
       subplot(3,1,1);
       imagesc(rho(2:numofgrids-1,2:numofgrids-1))
%        surf(X,Y,rho)
%        view(2)
%        xlabel('mm')
%        ylabel('mm')
%        zlabel('rho: active receptor ratio')
%        title(['Time frame #: ', num2str(i)])
%        colorbar
       
       subplot(3,1,2);
       imagesc(beta(2:numofgrids-1,2:numofgrids-1))
%        surf(X,Y,beta)
%        view(2)
%        xlabel('mm')
%        ylabel('mm')
%        zlabel('beta: [cAMP]in/K_R')
%        title(['Time frame #: ', num2str(i)])
%        colorbar
       
       
       subplot(3,1,3);
       imagesc(gamma(2:numofgrids-1,2:numofgrids-1))
%        surf(X,Y,gamma)
%        view(2)
%        xlabel('mm')
%        ylabel('mm')
%        zlabel('gamma: [cAMP]ex/K_R')
%        title(['Time frame #: ', num2str(i)])
%        colorbar
       
%             imagesc(gamma(2:numofgrids-1,2:numofgrids-1))
%             annotation(figure(2),'ellipse',...
%                 [0.282 0.82 0.02 0.02],'FaceColor',[1 0 0],...
%                 'Color',[1 0 0]);
      
       
%             imagesc(rho(2:numofgrids-1,2:numofgrids-1))
%             annotation(figure(2),'ellipse',...
%                 [0.282 0.348 0.02 0.02],'FaceColor',[1 0 0],...
%                 'Color',[1 0 0]);
%             
       set(f2, 'Position', [400 400 560 420])
       pause(0.06)
       
       % gamma and rho over time for a fixed point in space
       f3 = figure(3);
       plot(gamma(20,20),rho(20,20),'o') % set position of the observation point
      %  axis([-0.1 1.1 -0.1 1.1])
       hold off        
   end
   
   % Adding some extra v at time 25 and 50 for enjoyment
%    if ( 0 == mod(i,500) && add_extra )       
%        rho(45:55,45:55) = 0.9;       
%    end
   
   % Displaying elapsed time units
   clc; disp('Time: '); disp(t(i))      
end

% Adding labels and dots
figure(2)
subplot(2,1,1)
title('gamma [cAMP]ex/K_R','FontSize', 11, 'FontWeight', 'bold')
xlabel('x','FontWeight','bold','FontSize',14)
ylabel('y','FontWeight','bold','FontSize',14)
colorbar()

subplot(2,1,2)
title('rho: active receptor ratio','FontSize',11,'FontWeight','bold')
xlabel('x','FontWeight','bold','FontSize',14)
ylabel('y','FontWeight','bold','FontSize',14)
colorbar

figure(3)
xlabel('gamma','FontWeight','bold','FontSize',16)
ylabel('rho','FontWeight','bold','FontSize',16)

figure(4)
plot(2:numofgrids-1,gamma(2:numofgrids-1,20),'k',2:numofgrids-1,rho(2:numofgrids-1,20),'r')
xlabel('y','FontWeight','bold','FontSize',16)
ylabel('Concentration','FontWeight','bold','FontSize',16)
legend('gamma', 'rho')
