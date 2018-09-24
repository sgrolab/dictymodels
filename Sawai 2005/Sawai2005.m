%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is based on Saiwai 2005 model. The simulation is run on a
% square meshgrid, with each aquare represents a cell (if there is a cell).
% The inputs are initial vectorized E (excitability),C (extracellular cAMP
% concentration), C_T(excitation threshold),state (state of the cell on the
% grid).Parameters in the model are either declared as global variable or
% function inputs.
%
% The cells can have three states: state 0 (excitable), state 1 (excited,
% either releasing cAMP or in absolute refractory period), and state 2
% (relative refractory period). The function goes over all the grid squares
% in the field to update C, E, C_T, release_mask, and state with every
% simulation time step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,E,C_T,state,release_mask]=Sawai2005(eta,beta,dt,t,C,E,C_T,state,release_mask)

global my_dir firing_ratio D h alpha gamma C_max C_min x y X Y w l conc_release A E_min
global cell_mask E_max num_of_dt_release T_ARP T_RRP C0
filename=strcat('FR',num2str( firing_ratio),' dt',num2str(dt),'D',num2str(D),' Emin',num2str(E_min),' eta',num2str(eta),' beta',num2str(beta),'dt_C0=',num2str(C0));
figure

for i=2:1:length(t)%  go through all the time points length(t)
    E_prev=E(:,i-1);state_prev=state(:,i-1);release_mask_prev=release_mask(:,i-1);
    C_T_prev=C_T(:,i-1);C_prev=C(:,i-1);
    E_now=E(:,i);state_now=state(:,i);release_mask_now=release_mask(:,i);
    C_T_now=C_T(:,i);C_now=C(:,i);
    for j=1:1:length(cell_mask) % cycle through all the grids
        if cell_mask(j,1)==1 % where there is cell
             E_new=E_prev(j)+dt.*(eta+beta.*C_prev(j));
             if E_new>=E_max
                 E_now(j)=E_max;
             else
                 E_now(j)=E_new;
             end
            if state_prev(j)==0 % if previous state is 0 (excitable)
                prev_0=prev_time(state(j,1:(i-1)),0);% number of previous state 0 time steps
                if C_prev(j)>C_T_prev(j) || (prev_0*dt)>=15 % If C is greater than threshold or this cell have not been excited for 15 minutes
                    state_now(j)=1;
                    C_T_now(j)=C_max;release_mask_now(j)=0;
                else
                    state_now(j)=0;C_T_now(j)=C_min;release_mask_now(j)=0;
                end
  
            elseif state_prev(j)==1 % if previous state is 1 (excited state)
                prev_1=prev_time(state(j,1:(i-1)),1); % number of previous state 1 time steps
                if prev_1<= num_of_dt_release 
                    release_mask_now(j)=1;
                    state_now(j)=1;C_T_now(j)=C_max;
                elseif prev_1>num_of_dt_release  && prev_1<=(T_ARP+num_of_dt_release)
                    release_mask_now(j)=0;
                    state_now(j)=1;C_T_now(j)=C_max;
                else
                    state_now(j)=2;
                    C_T_now(j)=C_max;release_mask_now(j)=0;
                end
                
            elseif state_prev(j)==2
                chuchu=prev_time(state(j,1:(i-1)),2); % number of previous state 2 time steps
                if chuchu<T_RRP && C_prev(j)<C_T_prev(j)
                    C_T_now(j)=(C_max-A.*chuchu./(chuchu+T_ARP)).*(1-E_now(j));
                    state_now(j)=2;release_mask_now(j)=0;
                elseif chuchu<T_RRP && C_prev(j)>=C_T_prev(j)
                    state_now(j)=1;
                    C_T_now(j)=C_max;release_mask_now(j)=0;
                elseif chuchu>=T_RRP
                    state_now(j)=0;
                    C_T_now(j)=C_min; 
                    release_mask_now(j)=0;
                end
                
            end
        else % for grids without a cell
            E_now(j)=-1;
            C_T_now(j)=-1;
            state_now(j)=-1;
            release_mask_now(j)=0; 
        end
    end
   
   E(:,i)=E_now;state(:,i)=state_now;release_mask(:,i)=release_mask_now;
   C_T(:,i)= C_T_now;C(:,i)=C_now;
    
    C_mat=reshape(C_prev,[w,l]); % convert vector C into matrix form
    laplace_mat=4.*del2(C_mat,h); % laplacian
    laplace_vec=laplace_mat(:); % convert laplacian back into vector form
    C_now=C_prev+dt.*(D.*laplace_vec-gamma.*C_prev+conc_release.*release_mask_prev) ;
    C_now(C_now(:)<0)=0;% set all the negative value to zero
    % zero flux boundary condition with C_present matrix
    C_now=reshape(C_now,[w,l]); % convert vector C into matrix form
    C_now(1,:)=C_now(2,:);C_now(end,:)=C_now(end-1,:); 
    C_now(:,1)=C_now(:,2);C_now(:,end)=C_now(:,end-1);
    C_now=C_now(:); C(:,i)=C_now; % reshape C back into vector
    
    if ~isempty(find(C_prev>10000,1)) % for overshoot simulation
        disp(strcat('Overshoots for dt=',num2str(dt)))
        return
    end

    % plot in real time
    if (0 == mod(i,50)|| i==length(t) ) 
        C_plot=C_now;
        C_plot=reshape(C_plot,[w,l]);
        surf(x,y,C_plot)
        shading interp;  view(2);
        xlabel('mm'); ylabel('mm'); zlabel('extracellular cAMP')
        colorbar
        title(strcat(filename,' at frame #: ', num2str(i)))
        pause(0.0000001)
        % save images
        if i==length(t)
             saveas(gcf,strcat(my_dir,'outputs\',filename,' end.jpg'));
        end
     end
end
% save (strcat(my_dir,filename,'result.mat'), 'E','C' ,'C_T', 'state', 'release_mask') ;
disp(strcat('simulation ended for dt=',num2str(dt),', eta ',num2str(eta),' beta=',num2str(beta)))
close all

% save simulation result as video
frame_interval=10;
C_plot=C(:,1);
C_plot=reshape(C_plot,[w,l]);
surf(x,y,C_plot);
title('Extracellular cAMP at time point #: 1')
shading interp;  view(2);
zlabel('extracellular cAMP')
colorbar
I(1)=getframe(gcf);
for ii=frame_interval:frame_interval:length(t) 
    C_plot=C(:,ii);
    C_plot=reshape(C_plot,[w,l]);
    surf(x,y,C_plot);
    title(['Extracellular cAMP at time point #: ', num2str(ii)])
    shading interp;  view(2);
    zlabel('extracellular cAMP')
    colorbar
    I(ii./frame_interval)=getframe(gcf);
end
video = VideoWriter(strcat(my_dir,'outputs\',filename,'.avi')); %create the video object
open(video); %open the file for writing
writeVideo(video,I); %write the image to file
close(video); %close the file
end
