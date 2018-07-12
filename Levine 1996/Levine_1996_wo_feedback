%% Levine 1996 without feedback (constant E)
% Images and movies saved at a set path. Please change the path.
clear all; close all;
% Parameters
conc_release=300; t_release=1; gamma=8; 
E_max=0.93; E_min=0.1; % not mentioned in the paper?
alpha=0.0005; beta=1.24;
D=1; % diffusian coefficient
% time step and period, non-dimensional
dt=1;
t=0:dt:400;
num_of_dt_release=t_release./dt;
conc_release_dt=conc_release./(t_release./dt); % cAMP release every time step
T_ARP=8./dt; T_RRP=2./dt; % # of time steps in ARP and RRP
% T_ARP=2./dt; T_RRP=7./dt; % Sawai 2015
C_max=100; C_min=4; A=(T_RRP+T_ARP).*(C_max-C_min)./T_RRP;

% create square meshgrid and place cells
h=1;
[x,y]=meshgrid(0:h:150,0:h:150);
w=length(x(1,:));l=length(y(:,1));
X=x(:);Y=y(:); % reshape into a vector

% Initializations
C=1.*ones(length(Y),1); % extracellular cAMP concentration
cell_mask=zeros(length(Y),1);
% where there are cells
cell_ratio=0.98; % cell density
cell_num=floor(cell_ratio.*length(cell_mask));
cell_index=randperm(length(cell_mask), cell_num);
cell_mask(cell_index)=1;

E=E_min*ones(length(Y),1);E(cell_mask==0)=-1; % excitability, mark grids w/o cells as -1

C_T=C_min*ones(length(Y),1); C_T(cell_mask==0)=-1;% excitation threshold
release_mask=zeros(length(Y),1);release_mask(cell_mask==0)=-1; % initialize release to be 0
% Most cells  start at state 0, and a few cells start at state 1
state=zeros(length(Y),1);state(cell_mask==0)=-1; 
firing_ratio=0.001; firing_num=floor(firing_ratio.*sum(cell_mask));
index=find(cell_mask); select=index(randperm(length(index), firing_num)); state(select)=1;

figure % show the cell initial states
surf(x,y,reshape(state,[w,l]))
shading interp; view(2);
%% Start simulation
figure
for i=2:1:length(t)%  go through all the time points length(t)
    for j=1:1:length(cell_mask) % cycle through all the grids
        if cell_mask(j,1)==1 % where there is cell
            
%              E_new=E(j,i-1)+dt.*(-alpha.*E(j,i-1)+beta.*C(j,i-1));
%                 if E_new>E_max
%                     E(j,i)=E(j,i-1);
%                 else
%                     E(j,i)=E_new;
%                 end
                
            if state(j,i-1)==0 % previors state is 0
                prev_0=prev_time(state(j,1:(i-1)),0);% # of previous state 0
                if C(j,i-1)>C_T(j,i-1) || (prev_0*dt)>=15
                    state(j,i)=1;
                    C_T(j,i)=C_min;release_mask(j,i)=0;
                else
                    state(j,i)=0;C_T(j,i)=C_min;release_mask(j,i)=0;
                end
  
            elseif state(j,i-1)==1
                prev_1=prev_time(state(j,1:(i-1)),1);
                if prev_1<= num_of_dt_release 
                    release_mask(j,i)=1;
                    state(j,i)=1;C_T(j,i)=C_max;
                elseif prev_1>num_of_dt_release  && prev_1<=T_ARP
                    release_mask(j,i)=0;
                    state(j,i)=1;C_T(j,i)=C_max;
                else
                    state(j,i)=2;
                    C_T(j,i)=C_max;release_mask(j,i)=0;
                end
                
            elseif state(j,i-1)==2
                chuchu=prev_time(state(j,1:(i-1)),2);
                if chuchu<T_RRP && C(j,i-1)<C_T(j,i-1)
                    C_T(j,i)=(C_max-A.*chuchu./(chuchu+T_ARP)).*(1-E_max);
                    state(j,i)=2;release_mask(j,i)=0;
                elseif chuchu<T_RRP && C(j,i-1)>=C_T(j,i-1)
                    state(j,i)=1;
                    C_T(i,j)=C_max;release_mask(j,i)=1;
                elseif chuchu>=T_RRP
                    state(j,i)=0;
                    C_T(j,i)=C_min; % should it jump back to C_min?
                    release_mask(j,i)=0;
                end
                
            end
        else % for grids w/o cell
            % E(j,i)=-1;
            C_T(j,i)=-1;
            state(j,i)=-1;
            release_mask(j,i)=0; % where there is no cell, there is no release
        end
    end
    C_mat=reshape(C(:,i-1),[w,l]);
    laplace_mat=4.*del2(C_mat,h); % laplacian
    laplace_vec=laplace_mat(:);
    C_present=C(:,i-1)+dt.*(D.*laplace_vec-gamma.*C(:,i-1)+conc_release_dt.*release_mask(:,i-1)) ;
    C_present(C_present(:)<0)=0;% set all the negative value to zero
    % zero flux boundary condition with C_present matrix
    C_present=reshape(C_present,[w,l]);
    C_present(1,:)=C_present(2,:);C_present(end,:)=C_present(end-1,:); 
    C_present(:,1)=C_present(:,2);C_present(:,end)=C_present(:,end-1);
    % reshape back to vector
    C_present=C_present(:); C(:,i)=C_present; 

    % plot in real time
    if (0 == mod(i,10) ) 
        C_plot=C(:,i);
        C_plot=reshape(C_plot,[w,l]);
        surf(x,y,C_plot)
        title(['Extracellular cAMP at time point #: ', num2str(i)])
        shading interp;  view(2);
        xlabel('mm'); ylabel('mm'); zlabel('extracellular cAMP')
        colorbar
        pause(0.001)
    end
end

