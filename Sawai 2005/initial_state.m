function [cell_mask,release_mask,state,X,Y,w,l,x,y]=initial_state(h,firing_ratio)
    % create square meshgrid and place cells
    [x,y]=meshgrid(0:h:150*h,0:h:150*h);
    w=length(x(1,:));l=length(y(:,1));
    X=x(:);Y=y(:); % reshape into a vector

    % Initializations
    cell_mask=zeros(length(Y),1);
    % where there are cells
    cell_ratio=0.98; % cell density
    cell_num=floor(cell_ratio.*length(cell_mask));
    cell_index=randperm(length(cell_mask), cell_num);
    cell_mask(cell_index)=1;

    release_mask=zeros(length(Y),1);release_mask(cell_mask==0)=-1; % initialize release to be 0
    %  Most cells  start at state 0, and a few cells start at state 1
    state=zeros(length(Y),1);state(cell_mask==0)=-1; 
    firing_num=floor(firing_ratio.*sum(cell_mask));
    index=find(cell_mask); select=index(randperm(length(index), firing_num)); state(select)=1;

    % save initial cell states
    Initial_data_name=strcat('U:\cilse_research_sgro\Chuqiao\cAMP modeling\Sawai 2005\initial_cond\initial_FR_',num2str(firing_ratio),'.mat');
    save(Initial_data_name, 'cell_mask', 'release_mask', 'state','X','Y','x','y','w','l'); 

    figure % show the cell initial states
    surf(x,y,reshape(state,[w,l]))
    shading interp; view(2);
    initial_fig_name=strcat('U:\cilse_research_sgro\Chuqiao\cAMP modeling\Sawai 2005\initial_cond\initial_FR_',num2str(firing_ratio),'.jpg');
    saveas(gcf,initial_fig_name);
end
