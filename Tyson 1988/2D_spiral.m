%% 2D spiral wave generation

clear all, clc, close all

add_extra=true; % extra stimulation

% Some parameters based on Table 2, parameter set C
L1=10; L2=0.005; kappa=18.5; c=10;  
lambda1=1e-3;
lambda2=2.4; 
s1=950; 
s2=0.05; 
epsilon=0.01;
epsilon_prime=0.0182 ;  
D=0.024;% mm2/min
kt=3.0;
ke=12;
k1=0.12; 
k_1=1.2;
k2=2.22;
k_2=0.011;
h=5;
% parametr set A (Goldbeter oscillation)
% L1=10; L2=0.005; kappa=18.5; c=10;  
% lambda1=1e-4; % 1e-3;
% lambda2=0.26; %2.4; 
% s1=690;% 950; 
% s2=0.033;% 0.05; 
% epsilon_prime=0.014; % 0.0182 ;  
% D=0.024; % mm2/min
% kt=0.9; % 3.0;
% ke=5.4; % 12;
% k1=0.036; %0.12; 
% k_1=0.36; % 1.2;
% k2=0.666;% 2.22;
% k_2=0.0033;% 0.011;
% h=5;
% Initializations
[X,Y]=meshgrid(0:0.045:9,0:0.045:9); %% field of observation 9mm X 9mm
 X=k1/sqrt(ke*D)*X; Y=k1/sqrt(ke*D)*Y; % dimensionless X and Y
numofgrids=length(X(1,:));
step=X(1,2)-X(1,1);
% Time step size
dt = 0.12/60*k1; % min
% Time vector
t= 0:dt:90;

% Initializing gamma beta and rho
rho = 0.1.*ones(numofgrids,numofgrids);
beta = zeros(numofgrids,numofgrids);
gamma = zeros(numofgrids,numofgrids);

%% randomly assign high and zero intracellular cAMP in the middle square- does not generate spiral pattern
cAMP_high=10;
x=linspace(0,2,numofgrids);
y=cAMP_high./(1+exp(5.*(-x+1)));
beta=repmat(y,numofgrids,1)
% surf(X,Y,beta)


% plot initial conditions
f1 = figure(1)
       subplot(2,2,1)
       imagesc(rho(2:numofgrids-1,2:numofgrids-1))
       title('Initial condition of rho','FontSize', 11, 'FontWeight', 'bold')
       subplot(2,2,2)        
       imagesc(beta(2:numofgrids-1,2:numofgrids-1))
       colorbar
       title('Initial condition of beta','FontSize', 11, 'FontWeight', 'bold')
       subplot(2,2,3)        
       imagesc(gamma(2:numofgrids-1,2:numofgrids-1))
       title('Initial condition of gamma','FontSize', 11, 'FontWeight', 'bold')
movegui(f1,'northwest')

% Start simulation

f1gamma=@(gamma) (1+kappa.*gamma)./(1+gamma);
f2gamma=@(gamma) (L1+kappa.*L2.*c.*gamma)./(1+c.*gamma);
phi=@(gamma,rho) (lambda1+(rho.*gamma./(1+gamma)).^2)./(lambda2+(rho.*gamma./(1+gamma)).^2);

for i = 1:1:length(t)
   % zero-flux boundaries
   gamma(2,1:end)=gamma(1,1:end);  gamma(end-1,1:end)=gamma(end,1:end);
   gamma(1:end,2)=gamma(1:end,1);  gamma(1:end,end-1)=gamma(1:end,end);
   
     L=4*del2(gamma,step);
     
     rho_new= rho+dt.*(-f1gamma(gamma).*rho+f2gamma(gamma).*(ones(numofgrids,numofgrids)-rho));
     beta_new= beta+dt./epsilon_prime.*(s1.*phi(gamma,rho)-beta);
     gamma_new = gamma+dt.*((s2.*beta-gamma)./epsilon+epsilon.*L);
     
%      rho_new= rho+dt.*k1.*(-f1gamma(gamma).*rho+f2gamma(gamma).*(ones(numofgrids,numofgrids)-rho));
%      beta_new= beta+dt.*k1./epsilon_prime.*(s1.*phi(gamma,rho)-beta);
%      gamma_new = gamma+dt.*(kt./h.*beta-ke.*gamma+D.*L);
   rho=rho_new;
   beta=beta_new;
   gamma=gamma_new;
   
   % Plotting every tenth time step
   if (0 == mod(i,50) )       
       f2 = figure(2);
%        subplot(2,2,1)
%        surf(X,Y,rho)
%        shading interp; view(2);
%        xlabel('mm'); ylabel('mm'); zlabel('rho: active receptor ratio')
%        title(['Active receptor ratio at time point #: ', num2str(i)])
%        colorbar
       
%        subplot(2,2,2)
%        surf(X,Y,beta)
%        shading interp; view(2);
%        xlabel('mm'); ylabel('mm'); zlabel('beta: intracellular cAMP')
%        title(['Intracellular cAMP at time point #: ', num2str(i)])
%        colorbar
       
%        subplot(2,2,3)
       surf(X,Y,gamma)
       shading interp; view(2);
       xlabel('mm'); ylabel('mm'); zlabel('gamma: extracellular cAMP')
       title(['Extracellular cAMP at time point #: ', num2str(i)])
       colorbar
       
       pause(0.00000001)
   end
   
   % Adding some extra stimulation
%    if ( i<200 && add_extra )       
%        gamma(y1:(y2-2),x1:(x2-1))=10.*ones((l-1),(w));     
%    end
%    
   % Displaying elapsed time units
   clc; disp('Time: '); disp(t(i))      
end
