% Gregor 2010 model- single cell behavior

%% Dose response of oscillation period to external cAMP
clear all
clc
dt=0.01;
k=5;% ml/min
noise=0.002; % noise larger than in the paper

t1=0:dt:1000;
camp_cyt0=0.4;
% camp_ex_input=10e-6; % uM, oscillation threshold 5e-6~ 10e-6 uM
camp_ex_input=logspace(-5,-3,20);
for i=1:length(camp_ex_input)
    [camp_cyt1(i,:),camp_ex1,thetai1,noise1]=single_cell_fun(camp_cyt0,camp_ex_input(i),t1,noise);
    [PKS,LOCS,W]=findpeaks(camp_cyt1(i,:),t1,'WidthReference','halfheight','MinPeakHeight',0.5,'MinPeakDistance',0.5,'MinPeakProminence',0.5);
    Osctime=diff(LOCS);
    Osctime(isnan(Osctime))=0;Osctime(isempty(Osctime))=0;
    MeanOscPk(i)=mean(Osctime);
end

plot(camp_ex_input, MeanOscPk)
xlabel('cAMP_{ext}/uM');ylabel('Mean ocillation period/ min')
title('Dose response of oscillation period to cAMP_{ext}')

