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
% in the field to update C, E, C_T, release_mask, and statewiith every
% simulation time step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,E,C_T,state,release_mask]=Sawai2005(eta,beta,dt,t,C,E,C_T,state,release_mask)

global my_dir firing_ratio D h alpha gamma C_max C_min x y X Y w l conc_release A E_min
filename=strcat('FR',num2str( firing_ratio),' dt',num2str(dt),' Emin',num2str(E_min),' eta',num2str(eta),' beta',num2str(beta),'dt_C0=0.1');
global cell_mask E_max num_of_dt_release T_ARP T_RRP
figure

for i=2:1:length(t)%  go through all the time points length(t)
    for j=1:1:length(cell_mask) % cycle through all the grids
        if cell_mask(j,1)==1 % where there is cell
             E_new=E(j,i-1)+dt.*(eta+beta.*C(j,i-1));
             if E_new>=E_max
                 E(j,i)=E_max;
             else
                 E(j,i)=E_new;
             end
            if state(j,i-1)==0 % if previous state is 0 (excitable)
                prev_0=prev_time(state(j,1:(i-1)),0);% number of previous state 0 time steps
                if C(j,i-1)>C_T(j,i-1) || (prev_0*dt)>=15 % If C is greater than threshold or this cell have not been excited for 15 minutes
                    state(j,i)=1;
                    C_T(j,i)=C_max;release_mask(j,i)=0;
                else
                    state(j,i)=0;C_T(j,i)=C_min;release_mask(j,i)=0;
                end
  
            elseif state(j,i-1)==1 % if previous state is 1 (excited state)
                prev_1=prev_time(state(j,1:(i-1)),1); % number of previous state 1 time steps
                if prev_1<= num_of_dt_release 
                    release_mask(j,i)=1;
                    state(j,i)=1;C_T(j,i)=C_max;
                elseif prev_1>num_of_dt_release  && prev_1<=(T_ARP+num_of_dt_release)
                    release_mask(j,i)=0;
                    state(j,i)=1;C_T(j,i)=C_max;
                else
                    state(j,i)=2;
                    C_T(j,i)=C_max;release_mask(j,i)=0;
                end
                
            elseif state(j,i-1)==2
                chuchu=prev_time(state(j,1:(i-1)),2); % number of previous state 2 time steps
                if chuchu<T_RRP && C(j,i-1)<C_T(j,i-1)
                    C_T(j,i)=(C_max-A.*chuchu./(chuchu+T_ARP)).*(1-E(j,i));
                    state(j,i)=2;release_mask(j,i)=0;
                elseif chuchu<T_RRP && C(j,i-1)>=C_T(j,i-1)
                    state(j,i)=1;
                    C_T(i,j)=C_max;release_mask(j,i)=0;
                elseif chuchu>=T_RRP
                    state(j,i)=0;
                    C_T(j,i)=C_min; 
                    release_mask(j,i)=0;
                end
                
            end
        else % for grids without a cell
            E(j,i)=-1;
            C_T(j,i)=-1;
            state(j,i)=-1;
            release_mask(j,i)=0; 
        end
    end
   
    C_mat=reshape(C(:,i-1),[w,l]); % convert vector C into matrix form
    laplace_mat=4.*del2(C_mat,h); % laplacian
    laplace_vec=laplace_mat(:); % convert laplacian back into vector form
    C_present=C(:,i-1)+dt.*(D.*laplace_vec-gamma.*C(:,i-1)+conc_release.*release_mask(:,i-1)) ;
    C_present(C_present(:)<0)=0;% set all the negative value to zero
    % zero flux boundary condition with C_present matrix
    C_present=reshape(C_present,[w,l]); % convert vector C into matrix form
    C_present(1,:)=C_present(2,:);C_present(end,:)=C_present(end-1,:); 
    C_present(:,1)=C_present(:,2);C_present(:,end)=C_present(:,end-1);
    C_present=C_present(:); C(:,i)=C_present; % reshape C back into vector
    
    if ~isempty(find(C(:,i)>10000,1)) % for overshoot simulation
        disp(strcat('Overshoots for dt=',num2str(dt)))
        return
    end

    % plot in real time
    if (0 == mod(i,10)|| i==length(t) ) 
        C_plot=C(:,i);
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

% % save simulation result as video
% frame_interval=10;
% C_plot=C(:,1);
% C_plot=reshape(C_plot,[w,l]);
% surf(x,y,C_plot);
% title('Extracellular cAMP at time point #: 1')
% shading interp;  view(2);
% zlabel('extracellular cAMP')
% colorbar
% I(1)=getframe(gcf);
% for ii=frame_interval:frame_interval:length(t) 
%     C_plot=C(:,ii);
%     C_plot=reshape(C_plot,[w,l]);
%     surf(x,y,C_plot);
%     title(['Extracellular cAMP at time point #: ', num2str(ii)])
%     shading interp;  view(2);
%     zlabel('extracellular cAMP')
%     colorbar
%     I(ii./frame_interval)=getframe(gcf);
% end
% video = VideoWriter(strcat(my_dir,'outputs\',filename,'.avi')); %create the video object
% open(video); %open the file for writing
% writeVideo(video,I); %write the image to file
% close(video); %close the file
end
