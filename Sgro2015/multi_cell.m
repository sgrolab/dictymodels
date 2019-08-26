% multi_cell.m describes the population model in Sgro 2015. 
% INPUTS:
% alphaf: cAMP flow-in rate 
% j: cAMP out flow rate 
% rho: cell density
% para: parameters in the model
% dt: time step size
% t: simulation time vector
% OUTPUTS:
% Outputs are matrix representing variable at specific time for all the
% cells in a population. At a single time point it is represented as column vector.
% A: activator concentration. 
% R: repressor concentration
% cAMP: external cAMP concentration
function [A,R,cAMP] =multi_cell (alphaf,j,rho,para,dt,t,time_separation)
e=para.e; %tauA/tauR: episilon, ratio between activator and repressor time scales
tauA=para.tauA;
tauR=para.tauR;
g=para.g; %gamma, repressor degradation rate
c0=para.c0; % steady state repressor value in the absense of external cAMP
sigma=para.sigma; %sigma, noise strength
N=para.N; % number of cells in population simulations
a=para.a; % cAMP response magnitude
alpha0=para.alpha0;% basal cAMP leakage rate
alpha_pde=para.alpha_pde; % basal PDE leakage rate
Kd=para.Kd; % cAMP response threshold
S=para.S; % firing cell cAMP release rate                         

fA = @(A,R,cAMP) (A- (A.^3)/3 - R +a.*log(1 + cAMP./Kd) ) ;            
fR= @(A,R) e.*(A-g.*R + c0);
fcAMP= @(A,cAMP) alphaf + rho.*alpha0 + rho.*S.*sum(heaviside(A))./N -(alpha_pde.*rho +j).*cAMP;                             

n = length(t);
A = zeros(N,n);
R = zeros(N,n);
cAMP=zeros(1,n);

A(:,1) =-3.*rand; 
R(:,1) =-1.*rand;
cAMP(1,1) =0;

stream=RandStream('mt19937ar','Seed',4);
RandStream.setGlobalStream(stream);
r = sqrt(dt)*randn(N,n-1);

for i = 1:n-1
    
    A(:,i+1)     = A(:,i)+fA(A(:,i),R(:,i),cAMP(:,i)).*dt+sigma*r(:,i);
    R(:,i+1)     = R(:,i)+ fR(A(:,i),R(:,i)).*dt;
    if time_separation==0
        cAMP(:,i+1)  = cAMP(:,i) + fcAMP(A(:,i),cAMP(:,i)).*dt ;
    else
        cAMP(:,i+1) = (alphaf + rho.*alpha0 + rho.*S.*sum(heaviside(A(:,i)))./N)./(alpha_pde.*rho +j);
    end
end


end