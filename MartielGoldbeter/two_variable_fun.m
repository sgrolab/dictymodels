function dxdt=two_variable_fun(t,x,alpha,sigma,q,ki,kt,h,ke,k1,k2,L1,L2,c,theta,lambda,eps,beta)
% 3 variable model
rho_T=x(1);
gamma=x(2);

f1gamma=(k1+k2.*gamma)./(1+gamma);
f2gamma=(k1.*L1+k2.*L2.*c.*gamma)./(1+c.*gamma);
Y=rho_T.*gamma./(1+gamma);
phi=alpha.*(lambda.*theta+eps.*Y.^2)./(1+alpha.*theta+eps.*Y.^2.*(1+alpha));
beta=q.*sigma.*phi./(kt+ki);
q_prime=q.*kt./(h.*(ki+kt));

% differential equations
dxdt_1=-f1gamma.*rho_T+f2gamma.*(1-rho_T);
dxdt_2=q_prime.*sigma.*phi-ke.*gamma;

dxdt=[dxdt_1; dxdt_2];
end
