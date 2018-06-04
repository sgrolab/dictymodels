clear; clc; close all;

j=1;
rho=1e-4;
alphaf=1e1;

[t,A,R,C] =MultiCell_func(alphaf,j,rho);

for i=1:100 
    plot(t,A(i,:)); hold on
end

plot(t,mean(A,1),'g','LineWidth',3)


