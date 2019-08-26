function [camp_cyt,camp_ex,thetai,r]=single_cell_fun (camp_cyt0,camp_ex_input,t,eta)
Amax=20; 
Abas=0.4; % uM
w=2*pi/6; % min-1
Vc=1.1e-9; % ml
St=1.33; % cm2
Sc=1.3e-6; % cm2
K=0.0004; % uM, 400 pM
c_sec= 3.6; % min-1
c_excite=1.01; % min-1
n=length(t);

% initial conditions
camp_cyt=zeros(1,n);
camp_cyt(1)=camp_cyt0;

camp_ex=camp_ex_input*ones(1,n);

sinthetai0=(camp_cyt0*2-Amax-Abas)./(-Amax+Abas);% Initial sin theta_i calculated from cAMP_cyt_i

thetai=zeros(1,n); 
thetai(1)=asin(sinthetai0); % convert vector from asin to phase angle

dt=t(2)-t(1);
% noise 
stream=RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(stream);
r = sqrt(dt)*randn(1,n-1);

for i=1:(n-1)
 % Runge Kutta
    k1=dt*(w*(1-K./(camp_ex(i)+K)*c_excite.*sin(thetai(1,i))));
    thetai_temp=thetai(1,i)+k1/2;
    k2=dt*(w*(1-K./(camp_ex(i)+K)*c_excite.*sin(thetai_temp)));
    thetai_temp=thetai(1,i)+k2/2;
    k3=dt*(w*(1-K./(camp_ex(i)+K)*c_excite.*sin(thetai_temp)));
    thetai_temp=thetai(1,i)+k3/2;
    k4=dt*(w*(1-K./(camp_ex(i)+K)*c_excite.*sin(thetai_temp)));
    thetai(1,i+1)=thetai(1,i)+1/6*(k1+2*k2+2*k3+k4)+eta*r(1,i);
%   euler's method
%    thetai(:,i+1)=thetai(:,i)+h*(w*(1-K./(camp_ex(i)+K)*c_excite.*sin(thetai(:,i))))+eta*r(:,i);
    camp_cyt(1,i+1)=((-Amax+Abas)*sin(thetai(1,i+1))+(Amax+Abas))/2;
    
%     % plot in real time 
%     if (0==mod(i,100))
%         subplot(3,1,1)
%        % plot cyoplasmic cAMP in real time
%        scatter(t(i),camp_cyt(1,i),'.')
%        xlabel('time/min')
%        ylabel('cAMP_{cyt}/\muM')
%        hold on
%        
%        subplot(3,1,2)
%        % plot d theta/dt
%        theta_plot=linspace(-pi/2,pi/2,100);
%        dtheta_dt=w*(1-K./(camp_ex(i)+K)*c_excite*sin(theta_plot));
%        plot(theta_plot,dtheta_dt,'b')
%        xlabel('theta')
%        ylabel('d\theta/dt')
%        xlim([-pi/2 pi/2])
%        ylim([-0.2 2])
%        
%        subplot(3,1,3)
%        % plot theta in real time
%        sin_thetai=sin(thetai(i));theta_range=asin(sin_thetai);% theta value within -pi/2 and pi/2
%        scatter(theta_range,w*(1-K./(camp_ex(i)+K)*c_excite*sin(theta_range)))
%        xlabel('theta')
%        ylabel('d\theta/dt')
%        xlim([-pi/2 pi/2])
%        ylim([-0.2 2])
%        
%  
%     end
%      pause(0.0001)
end
end