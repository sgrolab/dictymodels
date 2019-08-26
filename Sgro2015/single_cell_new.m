function [A,R,t,lineVal,A_orig,R_orig]=single_cell_new(cAMP,t_tot)
e=0.1; %tauA/tauR; %episilon=0.1
tauA=0.09;
tauR=tauA/e;
g=0.5; %gamma=0.5
c0=1.2;
sigma=0.15; %sigma
N=100; % number of cells in population simulations
a=0.058;
alpha0=800;
alpha_pde= 1000;
Kd= 1e-5;
S=1e6;                          

Nt=27; % normalization factor of t
Na=3.5;  % normalization factor of A
offset_A=1.5;

dt=0.005;
t = 0:dt:t_tot;   
tTrig=2.5./dt.*Nt; % at this time step trigger happens

% % fig. 3 step increase vs ramp increase
% tStepOff=tTrig+Nt.*5./dt;% at this time step trigger stops
% tRampOn=tStepOff+Nt.*5./dt; % at this time step ramp trigger happens
% tRampOff=tRampOn+Nt.*5./dt;% at this time step ramp trigger stops
% % gf = (cAMP/0.000001)^(1/(Nt.*5./dt)-1);
% % k = 1:(Nt.*5./dt);cAMP_ramp = 0.000001*gf.^(k-1);
% ln_cAMP_ramp=zeros(tRampOff-tRampOn+1);
% ln_cAMP_ramp(end)=log(cAMP);
% ln_cAMP_ramp(1)=log(0.01);
% ln_cAMP_ramp=linspace(ln_cAMP_ramp(1),ln_cAMP_ramp(end),tRampOff-tRampOn+1);
% cAMP_ramp=exp(ln_cAMP_ramp);
% 
% t_ramp_plot=1:1:tRampOff-tRampOn+1;
% figure
% plot(t_ramp_plot,cAMP_ramp)

lineVal=(tTrig*dt)./Nt;% Actual time when the step triggering happens
% the iteration # at which the the cAmp is introduced

fA = @(t,pr,st) (pr(1)- pr(1)^3/3 - pr(2) +a*log(1 + st/Kd) ) ;            
fR= @(t,pr) e*(pr(1)-g*pr(2) + c0);
                                          
n = length(t);
A = zeros(1,n);
R = zeros(1,n);

A(1) = -1.5; 
R(1) = -0.5;% Initial conditions read from fig 1E

% rng('default')
r = sqrt(dt)*randn(1,n-1)';

% % continuous step cAMP triggering
% for i = 1:n-1
%     if(i<tTrig)
%     A(i+1) = A(i)+fA(t(i),[A(i) R(i)],0)*dt+ sigma*r(i);
%     else
%     A(i+1) = A(i)+fA(t(i),[A(i) R(i)],cAMP)*dt+ sigma*r(i);
%     end
%     R(i+1) = R(i)+ fR(t(i),[A(i) R(i)])*dt/tauR;
% end

% constant followed by ramp cAMP triggering 
% for i = 1:n-1
%     if(i<tTrig)
%     A(i+1) = A(i)+fA(t(i),[A(i) R(i)],0)*dt+ sigma*r(i);
%     elseif(i>=tTrig && i<tStepOff)
%     A(i+1) = A(i)+fA(t(i),[A(i) R(i)],cAMP)*dt+ sigma*r(i);
%     elseif (i>=tStepOff && i<tRampOn)
%     A(i+1) = A(i)+fA(t(i),[A(i) R(i)],0)*dt+ sigma*r(i);
%     elseif (i>=tRampOn && i<tRampOff)
%     A(i+1) = A(i)+fA(t(i),[A(i) R(i)],cAMP_ramp(i-tRampOn+1))*dt+ sigma*r(i);
%     else
%     A(i+1) = A(i)+fA(t(i),[A(i) R(i)],0)*dt+ sigma*r(i);    
%     end
%     R(i+1) = R(i)+ fR(t(i),[A(i) R(i)])*dt/tauR;
% end
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
R=R./Na;
t=t./Nt;


end



