function [A,R,t,lineVal,A_orig,R_orig]=single_cell_old(cAMP,t_tot)
e=0.1; %tauA/tauR; %episilon
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

Nt=27; % normalization factor of A
Na=3.5;  % normalization factor of t
offset_A=1.5;
tTrig=5400;

dt=0.005;
t = 0:dt:t_tot;   

lineVal=(tTrig*dt)./Nt;
% the iteration # at which the the cAmp is introduced

fA = @(t,pr,st) (pr(1)- pr(1)^3/3 - pr(2) +a*log(1 + st/Kd) ) ;            
fR= @(t,pr) e*(pr(1)-g*pr(2) + c0);

                                            

n = length(t);
A = zeros(1,n);
R = zeros(1,n);

A(1) =- 1.5; 
R(1) = -0.5;

rng('default')
r = sqrt(dt)*randn(1,n-1)';

for i = 1:n-1
    if(i<tTrig)
    A(i+1) = A(i)+fA(t(i),[A(i) R(i)],0)*dt+ sigma*r(i);
    else
    A(i+1) = A(i)+fA(t(i),[A(i) R(i)],cAMP)*dt+ sigma*r(i);
    end
    R(i+1) = R(i)+ fR(t(i),[A(i) R(i)])*dt/tauR;
end
A_orig=A;R_orig=R;% un-scaled, un-offset A and R after simulation

A=(A+offset_A)./Na;
% R=R./Na;
t=t./Nt;



