% This function sets up initial states of the Levine model and save the
% initial state data as  .m file in a set directory. 
% INPUTS: 
% h- rectangular grid size
% cell_ratio- cell density
% firing_ratio- the ratio of cells that fire at time 1
% InitialDir- directory to save initial states
% OUTPUTS are vectors collapsed from square matrix:
% cell_mask- binary vector showing whether there is cell on a point
% release_mask- all zero vector, no cells release cAMP at time 1
% state- vector showing what state each cell starts at. State 0 means
% resting, state 1 means during firing or in absolute refractory period, state 2
% means in relative refractory period
% X- vector for all x-coordinates
% Y- vector for all y-coordinates
% x- meshgrid matrix for x coordinates
% y- meshgrid matrix for y coordinates
% w- length of X vetor
% l- length of Y vector
% initial_data_name- name of .m file that saves initial states
function [cell_mask,release_mask,state,X,Y,w,l,x,y,Initial_data_name]=initial_state(h,cell_ratio,firing_ratio,InitialDir)
    % create square meshgrid and place cells
    [x,y]=meshgrid(0:h:150*h,0:h:150*h);
    w=length(x(1,:));l=length(y(:,1));
    X=x(:);Y=y(:); % reshape into a vector

    % Initializations
    cell_mask=zeros(length(Y),1);
    % where there are cells
    cell_num=floor(cell_ratio.*length(cell_mask));
    cell_index=randperm(length(cell_mask), cell_num);
    cell_mask(cell_index)=1;

    release_mask=zeros(length(Y),1);release_mask(cell_mask==0)=-1; % initialize release to be 0
    %  Most cells  start at state 0, and a few cells start at state 1
    state=zeros(length(Y),1);state(cell_mask==0)=-1; 
    firing_num=floor(firing_ratio.*sum(cell_mask));
    index=find(cell_mask); select=index(randperm(length(index), firing_num)); state(select)=1;

    % save initial cell states in assigned directory
    Initial_data_name=strcat(InitialDir,'\InitialState_',num2str(firing_ratio),'.mat');
    save(Initial_data_name, 'cell_mask', 'release_mask', 'state','X','Y','x','y','w','l','h'); 

    figure % show the cell initial states
    surf(x,y,reshape(state,[w,l]))
    shading interp; view(2);
    initial_fig_name=strcat(InitialDir,'\InitialState_',num2str(firing_ratio),'.jpg');
    saveas(gcf,initial_fig_name);
end
