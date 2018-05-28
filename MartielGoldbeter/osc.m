function dydt = osc(t,y)

%parameters
k1 = 0.036;     % per min
k2 = 0.666;     % per min
L1 = 10;
L2 = 0.005 ;
c = 10;           % 0.15 - 50
lamda=0.01;
theta=0.01;
e=1;
q=4000;
sig=0.6;
v=12;
k= 4;         %k prime in the paper
ki=1.7 ;
kt=0.9;
kc=5.4;
h=5;


%rho = y(1);alpha=y(2);beta=y(3);gamma=y(4);

f1 =@(pr) (k1 +k2.*pr)./(1 + pr);
f2 =@(pr) (k1*L1 + k2*L2*c.*pr)./(1 + c.*pr);
Ysq  =@(pr1,pr2) (pr1.*pr2./(1+ pr2)).^2;
PI =@(pr1,pr2,pr3) pr3.*(lamda*theta + e.*Ysq(pr1,pr2))./(1 + theta.*pr3 + (1 + pr3).*e*Ysq(pr1,pr2) );

%dydt(1) is dy1/dt and etc..
dydt = zeros(4,1);    % this creates an empty column     %vector that you can fill with your four derivatives:    

dydt(1) = -f1(y(4)).*y(1) + f2(y(4)).*(1-y(1));    
dydt(2) = v - k*y(2) - sig*PI(y(1),y(4),y(2)); 
dydt(3) = q*sig.*PI(y(1),y(4),y(2)) - (ki + kt).*y(3);   
dydt(4) = (kt/h).*y(3) - kc.*y(4);   
end

