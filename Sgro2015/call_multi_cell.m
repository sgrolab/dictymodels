clear; clc; close all;

j=0:0.2:1;
logrho=-5.5:0.5:-3; %log of 10
alphafval=10;
% [x,y,z]=meshgrid(j,rho,alphaf);
n=length(j);
m=length(logrho);
%alphaf=alphafval.*ones(n,m);
rate=zeros(n,m);
itr=0;
for p=1:n
    for q=1:m
rate(q,p)=multi_cell (alphafval,j(p),10^(logrho(q)),0);
fprintf('p= %d of %d and q = %d of %d\n', p,n,q,m);
J(q,p)=j(p);
Rho(q,p)=logrho(q);

    end
end

map = jet(200);
rate=0.6.*rate./max(max(rate));
resizedRate = imresize(rate, 4, 'bicubic');
imshow(resizedRate,'colormap',map,'XData',j, 'YData', logrho); hold on
set(gca, 'ydir', 'normal')
colorbar
caxis([0 0.6])

axis square
axis on 