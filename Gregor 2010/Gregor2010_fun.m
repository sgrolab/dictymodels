 function [camp_cyt,camp_ex,thetai,r]=Gregor2010_fun (camp_cyt0,rho,k,t,eta)
Amax=20; 
Abas=0.4; % uM
w=2*pi/6; % min-1
Vc=1.1e-9; % ml
St=1.33; % cm2
Sc=1.3e-6; % cm2
K=0.0004; % uM, 400 pM
c_sec= 3.6; % min-1
c_excite=1.01; % min-1
Nc=length(camp_cyt0); % number of cells
n=length(t);

% initial conditions
camp_cyt=zeros(Nc,n);
camp_cyt(:,1)=camp_cyt0;

camp_ex=zeros(1,n);
camp_ex(1)=Vc*St/Sc*rho/k*c_sec*1/Nc*sum(camp_cyt(:,1));

sinthetai0=(camp_cyt0*2-Amax-Abas)./(-Amax+Abas);% Initial sin theta_i calculated from cAMP_cyt_i

thetai=zeros(Nc,n); 
thetai(:,1)=asin(sinthetai0); % convert vector from asin to phase angle

dt=t(2)-t(1);
% noise 
stream=RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(stream);
r = sqrt(dt)*randn(Nc,n-1);

for i=1:(n-1)
 % Runge Kutta
    k1=dt*(w*(ones(Nc,1)-K*ones(Nc,1)./(camp_ex(i)+K)*c_excite.*sin(thetai(:,i))));
    thetai_temp=thetai(:,i)+k1/2;
    k2=dt*(w*(ones(Nc,1)-K*ones(Nc,1)./(camp_ex(i)+K)*c_excite.*sin(thetai_temp)));
    thetai_temp=thetai(:,i)+k2/2;
    k3=dt*(w*(ones(Nc,1)-K*ones(Nc,1)./(camp_ex(i)+K)*c_excite.*sin(thetai_temp)));
    thetai_temp=thetai(:,i)+k3/2;
    k4=dt*(w*(ones(Nc,1)-K*ones(Nc,1)./(camp_ex(i)+K)*c_excite.*sin(thetai_temp)));
    thetai(:,i+1)=thetai(:,i)+1/6*(k1+2*k2+2*k3+k4)+eta*r(:,i);
%   euler's method
%    thetai(:,i+1)=thetai(:,i)+h*(w*(1-K./(camp_ex(i)+K)*c_excite.*sin(thetai(:,i))))+eta*r(:,i);
    camp_cyt(:,i+1)=((-Amax+Abas)*sin(thetai(:,i+1))+(Amax+Abas)*ones(Nc,1))/2;
    camp_ex(i+1)=Vc*St/Sc*rho/k*c_sec*1/Nc*sum(camp_cyt(:,i+1));
%     if (0==mod(i,100))
%        % plot 1st cell in real time
%        scatter(t(i),camp_cyt(1,i),'.')
%        hold on
%        pause(0.000001)
%     end  
end
end