

clc; clear all; close all
e=0.1; %tauA/tauR; %epicylon
tauA=0.2;
tauR=tauA/e;
g=0.58; %gamma
c0=1.2;
sigma=0.15; %sigma
N=100; % number of cells in population simulations
a=0.058;
alpha0=800;
alpha_pde= 1000;
Kd= 1e-5;
S=1e6;
cAmp=10;


fA = @(t,pr,st) (pr(1)- pr(1)^3/3 - pr(2) +a*log(1 + st/Kd) )/tauA  ;            
fR= @(t,pr) (pr(1)-g*pr(2) + c0)/tauR;

dt=0.005;
t = 0:dt:100;                   % Time vector
                             

n = length(t);
A = zeros(1,n);
R = zeros(1,n);

A(1) =- 1.5; 
R(1) = -0.5;



rng('default')
r = sqrt(dt)*randn(1,n-1)';


for i = 1:n-1
    if(i<4000)
    A(i+1) = A(i)+fA(t(i),[A(i) R(i)],0)*dt+ sigma*r(i);
    else
    A(i+1) = A(i)+fA(t(i),[A(i) R(i)],cAmp)*dt+ sigma*r(i);
    end
    R(i+1) = R(i)+ fR(t(i),[A(i) R(i)])*dt/tauR;
end

subplot (121); plot(t,A)
ylabel('A')
xlabel('t')
subplot (122); plot(t,R)
ylabel('R')
xlabel('t')

figure(2)
plot(A,R)
