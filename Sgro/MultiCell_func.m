function [t,A,R,Camp] =multi_cell (alphaf,j,rho)
e=0.1; %tauA/tauR; %epicylon
tauA=0.09;
tauR=tauA/e;
g=0.5; %gamma
c0=1.2;
sigma=0.15; %sigma
N=100; % number of cells in population simulations
a=0.058;
alpha0=800;
alpha_pde= 1000;
Kd= 1e-5;
S=1e6;                          



fA = @(A,R,Camp) (A- (A.^3)/3 - R +a.*log(1 + Camp./Kd) ) ;            
fR= @(A,R) e.*(A-g.*R + c0);
fCamp= @(A,Camp) alphaf + rho*alpha0 + rho*S.*sum(heaviside(A)) -(alpha_pde*rho +j).*Camp;


dt=0.005;
t = 0:dt:500;                   % Time vector
                             

n = length(t);
A = zeros(N,n);
R = zeros(N,n);

A(:,1) =-3.*rand; 
R(:,1) =-1.*rand;
Camp(:,1) =0;



rng(1,'twister')
r = sqrt(dt)*randn(N,n-1);


for i = 1:n-1
    
    A(:,i+1)     = A(:,i)+fA(A(:,i),R(:,i),Camp(:,i)).*dt+ sigma*r(:,i);
    R(:,i+1)     = R(:,i)+ fR(A(:,i),R(:,i)).*dt;
    Camp(:,i+1)  = Camp(:,i) + fCamp(A(:,i),Camp(:,i)).*dt ;
end


end