l = linspace(-10,10,100);
g = 30; 
y = l.^3+(5/3+g)*l.^2+ (2/3-5/51*g)*l + 38/51*g; 
xaxis = zeros(length(l));
figure()
plot(l,y)
hold on
plot(l,xaxis)
xline(0)
title(strcat('gamma = ',num2str(g)))
%ylim([-50,50])

%% 
g = 1000;
a = 1;
b = g+5/3;
c = 2/3-5/3*g;
d = 8/51*g;
x = roots([a b c d])