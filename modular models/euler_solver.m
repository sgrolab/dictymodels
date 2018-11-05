function [tout,xout]=euler_solver(fun,t,x0,param,stim)
    tout=t;
    odesize=length(x0); % order of system of odes
    % Check inputs
    if nargin < 5
      stim=zeros(size(t)); % zero extracellular cAMP stimulus at all time points
      if nargin < 4
            error(message('euler_colver:NotEnoughInputs'));
      end
    end
    n=length(t); dt=t(2)-t(1);
    xout(:,1)=x0;
    
    if param.sigma==[]
       sigma=zeros(size(x0));% no noise case
    end
    rng('default')
    r = sqrt(dt)*randn(1,n-1);
    r=[r;repmat(zeros(size(r)),odesize-1,1)]; % noise term added to the first ode equation, others set as zero
    for i=1:n-1
        dxdt=fun(xout(:,i),param,stim(i));
        xout(:,i+1)=xout(:,i)+dxdt.*dt+param.sigma*r(:,i);  
        % disp(strcat('A_next',num2str(xout(1,i+1))))
    end
end
