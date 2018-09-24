function [C,E,C_T,state,release_mask]=Levine1996(feedback,E_min,dt,t,C,E,C_T,state,release_mask)
global my_dir firing_ratio D h alpha beta gamma C_max C_min x y X Y w l conc_release A
filename=strcat('FR',num2str( firing_ratio),' dt',num2str(dt),' Emin',num2str(E_min),' FB',num2str(feedback));
global dE cell_mask E_max num_of_dt_release T_ARP T_RRP
figure
for i=2:1:length(t)%  go through all the time points length(t)
    for j=1:1:length(cell_mask) % cycle through all the grids
        if cell_mask(j,1)==1 % where there is cell
             if feedback==1 % with feedback
                 E_new=E(j,i-1)+dt.*(-alpha.*E(j,i-1)+beta.*C(j,i-1));
                    
             elseif feedback==2 % without feedback
                 E(j,i)=E_max; 
             else % with uniform feedack
                 E(j,i)=E(j,i-1)+dE;
             end
             if E_new>E_max
                    E(j,i)=E(j,i-1);
                else
                    E(j,i)=E_new;
             end
            if state(j,i-1)==0 % previors state is 0
                prev_0=prev_time(state(j,1:(i-1)),0);% # of previous state 0
                if C(j,i-1)>C_T(j,i-1) || (prev_0*dt)>=15
                    state(j,i)=1;
                    C_T(j, i)=C_min;release_mask(j,i)=0;
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
                    C_T(j,i)=(C_max-A.*chuchu./(chuchu+T_ARP)).*(1-E(j,i));
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
            E(j,i)=-1;
            C_T(j,i)=-1;
            state(j,i)=-1;
            release_mask(j,i)=0; % where there is no cell, there is no release
        end
    end
    C_mat=reshape(C(:,i-1),[w,l]);
    laplace_mat=4.*del2(C_mat,h); % laplacian
    laplace_vec=laplace_mat(:);
    C_present=C(:,i-1)+dt.*(D.*laplace_vec-gamma.*C(:,i-1)+conc_release.*release_mask(:,i-1)) ;
    C_present(C_present(:)<0)=0;% set all the negative value to zero
    % zero flux boundary condition with C_present matrix
    C_present=reshape(C_present,[w,l]);
    C_present(1,:)=C_present(2,:);C_present(end,:)=C_present(end-1,:); 
    C_present(:,1)=C_present(:,2);C_present(:,end)=C_present(:,end-1);
    % reshape back to vector
    C_present=C_present(:); C(:,i)=C_present; 
    
    if ~isempty(find(C(:,i)>10000)) % for overshoot simulation
        disp(strcat('Overshoots for dt=',num2str(dt),', feedback= ',num2str(feedback)))
        return
    end

    % plot in real time
    if (0 == mod(i,10)|| i==length(t) ) 
        C_plot=C(:,i);
        C_plot=reshape(C_plot,[w,l]);
        surf(x,y,C_plot)
        title(strcat(filename,' at frame #: ', num2str(i)))
        shading interp;  view(2);
        xlabel('mm'); ylabel('mm'); zlabel('extracellular cAMP')
        colorbar
        pause(0.0001)
        % save images
        if i==length(t)
             saveas(gcf,strcat(my_dir,filename,' end.jpg'));
        end
    end
end
save (strcat(my_dir,filename,'result.mat'), 'E','C' ,'C_T', 'state', 'release_mask') ;
disp(strcat('simulation ended for dt=',num2str(dt),', feedback= ',num2str(feedback)))
close all

% save as video
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
video = VideoWriter(strcat(my_dir,filename,'.avi')); %create the video object
open(video); %open the file for writing
writeVideo(video,I); %write the image to file
close(video); %close the file
end