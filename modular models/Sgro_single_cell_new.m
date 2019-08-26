function dxdt=Sgro_single_cell_new(x,param,stim)

a0=param.a0;
e=param.e; %tauA/tauR; %episilon
tauA=param.tauA;
tauR=param.tauR;
g=param.g; %gamma
c0=param.c0;
a=param.a;
alpha0=param.alpha0;
alpha_pde= param.alpha_pde;
Kd=param.Kd;
S=param.S;                          

dxdt(1,1) =(x(1)- x(1)^3/3 - 1/a0*x(2) +a*log(1 + stim/Kd) ) ;

% disp(strcat('dAdt=',num2str(dxdt(1))))
dxdt(2,1)= (x(1)-g*x(2) + c0).*e./tauR;

end
