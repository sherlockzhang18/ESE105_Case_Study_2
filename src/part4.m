clc;
close all;
clearvars;
% *******************************************************************
%% load the data
load("mockdata2023.mat");
plot([cumulativeDeaths;newInfections].');
% *******************************************************************
%% Before roll out
% initialize the initial state
x0 = [1-newInfections(1)-cumulativeDeaths(1),newInfections(1),0,0,cumulativeDeaths(1),0];
% parameters for the matrix
v_rate = 0;
inf_rate1 = 0.005;
inf_rate2 = 0.2;
immune_rate = 0.05;
immune_rate2 = 0.1;
death1 = 0.02;
death2 = 0.00005;
recover1 = 0.1;
recover2 = 0.07;
% update matrix
A = [1-inf_rate1-v_rate recover1                      0                              0            0 0
     inf_rate1          1-recover1-death1-immune_rate 0                              0            0 0
     0                  0                             1-recover2-death2-immune_rate2 0            0 inf_rate2
     0                  immune_rate                   immune_rate2                   1            0 0
     0                  death1                        death2                         0            1 0
     v_rate             0                             recover2                       0            0 1-inf_rate2 ];
% initialize time span
B = zeros(6,1);
t = 1:120;
% use the function to model
y = modelVaccine(A,B,t,x0);
% calculate the new cases from the model
model_newinf = ([x0(2)+x0(3),(y(:,1)*inf_rate1 + y(:,6)*inf_rate2).']);
model_newinf = model_newinf(1:120);
% *******************************************************************
%% After roll out
% change the values for the matrix
v_rate = 0.045;
inf_rate1 = 0.005;
inf_rate2 = 0.15;
immune_rate = 0.05;
immune_rate2 = 0.1;
death1 = 0.02;
death2 = 0.00005;
recover1 = 0.1;
recover2 = 0.07;

A = [1-inf_rate1-v_rate recover1                      0                              0           0 0
     inf_rate1          1-recover1-death1-immune_rate 0                              0           0 0
     0                  0                             1-recover2-death2-immune_rate2 0           0 inf_rate2
     0                  immune_rate                   immune_rate2                   1           0 0
     0                  death1                        death2                         0           1 0
     v_rate             0                             recover2                       0           0 1-inf_rate2 ];
% update the inital state and the time span
x0 = A*y(120,:).';
B = zeros(6,1);
t = (121:400)-120;
% call the function and model.
y1 = modelVaccine(A,B,t,x0);
model_newinf2 = (y1(:,1)*inf_rate1 + y1(:,6)*inf_rate2).';
% *******************************************************************
%% Combine the model and give output
% combine the model and the model-based newinfection.
x = [y;y1];
model_newinfection = [model_newinf,model_newinf2];
% use immse function to callibrate the parameters.
err1 = immse(x(:,5).',cumulativeDeaths);
err2 = immse(model_newinfection,newInfections);
% plot the figure 
figure;
plot([x(:,5).';model_newinfection].','LineWidth',2);
hold on;
plot([cumulativeDeaths;newInfections].','LineWidth',2)
title('Model with Vaccine')
legend('cumulative deaths (model)','new infection (model)','cumulative deaths','new infection','Location','northwest');
% calculate the vaxpop and vaxbreak and export
vaxpop = cumsum([y(:,1)*0;y1(:,1)*v_rate]);
vaxbreak = x(:,3);
save("competiton.mat",'vaxpop','vaxbreak');
% *******************************************************************
%% Define functions
function Y = modelVaccine(A,B,t,x0)
size1 = size(t,2);
sys_sir_base = ss(A,B,eye(6),zeros(6,1),1);
Y = lsim(sys_sir_base,zeros(size1,1),t,x0);
end