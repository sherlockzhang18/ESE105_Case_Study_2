clc;
clearvars;
close all;
% *******************************************************************
%% load and modify the data
load("COVID_STL.mat");
measured_deaths = deaths_STL/POP_STL;
measured_infected = cases_STL/POP_STL;
plot(measured_infected.');

plot(measured_deaths.');
% *******************************************************************
%% Pre-Delta Period
% From '2020-03-18' to '2021-06-23'

% set the intial state
x0 = [1-7/POP_STL, 7/POP_STL, 0, 0];
% matrix for the function
B = zeros(4,1);

%initialize the A matrix to update the system.
infectious_rate = 0.0015;
immune_rate = 0.39;
death_rate = 0.016;
recover_rate = 0.3;
rein_rate = 0.005;

A = [1-infectious_rate  recover_rate                            rein_rate   0;
     infectious_rate    1-(recover_rate+immune_rate+death_rate) 0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];

% identify the time span
t = 1:67;
% call the funtion and simulate the model
y = model(A,B,t,x0);
% calculate the model based cumulative cases per day.
new_cases_model = cumsum([x0(2),(y(:,1)*infectious_rate).']);
new_cases_model = new_cases_model(1:67);
% use immse function to callibrate the parameters.
err1 = immse(y(:,4).',measured_deaths(1:67));
err2 = immse(new_cases_model,measured_infected(1:67));

% plot the figures.
figure
plot(new_cases_model);
hold on;
plot(measured_infected(1:67));
plot(y(:,4));
plot(measured_deaths(1:67));
legend('cases (model)','measured infected','deceased (model)','measured deceased','location','northwest')
title('Pre-Delta Period')
% *******************************************************************
%% Delta-Period
% From '2021-06-30' to '2021-10-26'

%initialize the A matrix to update the system.
infectious_rate = 0.0017;
immune_rate = 0.35;
death_rate = 0.0075;
recover_rate = 0.35;
rein_rate = 0.005;

A = [1-infectious_rate  recover_rate                            rein_rate   0;
     infectious_rate    1-(recover_rate+immune_rate+death_rate) 0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];
% identify the time span
t = (68:84)-67;
% set the intial state
x0 = A*y(size(y,1),:).';
% call the funtion and simulate the model
y1 = model(A,B,t,x0);
% calculate the model based cumulative cases per day.
new_cases_model2 = new_cases_model(length(new_cases_model))+cumsum((y1(:,1)*infectious_rate).');
% use immse function to callibrate the parameters.
err3 = immse(y1(:,4).',measured_deaths(68:84));
err4 = immse(new_cases_model2,measured_infected(68:84));
% plot the figures.
figure
plot(new_cases_model2);
hold on;
plot(measured_infected(68:84));
plot(y1(:,4));
plot(measured_deaths(68:84));
legend('cases (model)','measured infected','deceased (model)','measured deceased','location','northwest')
title('Delta Period')
% *******************************************************************
%% Omicron Period
% From '2021-10-27' to '2022-03-22'

%initialize the A matrix to update the system.
infectious_rate = 0.0044;
immune_rate = 0.25;
death_rate = 0.005;
recover_rate = 0.4;
rein_rate = 0.05;
A = [1-infectious_rate  recover_rate                            rein_rate   0;
     infectious_rate    1-(recover_rate+immune_rate+death_rate) 0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];

% identify the time span
t = (85:105)-84;
% set the intial state
x0 = A*y1(size(y1,1),:).';
% call the funtion and simulate the model
y2 = model(A,B,t,x0);
% calculate the model based cumulative cases per day.
new_cases_model3 = new_cases_model2(length(new_cases_model2))+cumsum((y2(:,1)*infectious_rate).');
% use immse function to callibrate the parameters.
err5 = immse(y2(:,4).',measured_deaths(85:105));
err6 = immse(new_cases_model3,measured_infected(85:105));
% plot the figures
figure
plot(new_cases_model3);
hold on;
plot(measured_infected(85:105));
plot(y2(:,4));
plot(measured_deaths(85:105));
legend('cases (model)','measured infected','deceased (model)','measured deceased','location','northwest')
title('Omicron Period')
% *******************************************************************
%% Post-Omicron Period
% From '2021-10-27' to '2022-03-22'

%initialize the A matrix to update the system.
infectious_rate = 0.0017;
immune_rate = 0.2;
death_rate = 0.0045;
recover_rate = 0.5;
rein_rate = 0.05;
A = [1-infectious_rate  recover_rate                            rein_rate   0;
     infectious_rate    1-(recover_rate+immune_rate+death_rate) 0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];
% identify the time span
t = (106:158)-105;
% set the intial state
x0 = A*y2(size(y2,1),:).';
% call the funtion and simulate the model
y3 = model(A,B,t,x0);
% calculate the model based cumulative cases per day.
new_cases_model4 = new_cases_model3(length(new_cases_model3))+cumsum((y3(:,1)*infectious_rate).');
% use immse function to callibrate the parameters.
err7 = immse(y3(:,4).',measured_deaths(106:158));
err8 = immse(new_cases_model4,measured_infected(106:158));
% plot the figures.
figure
plot(new_cases_model4);
hold on;
plot(measured_infected(106:158));
plot(y3(:,4));
plot(measured_deaths(106:158));
legend('cases (model)','measured infected','deceased (model)','measured deceased','location','northwest')
title('Post-Omicron Period')
% *******************************************************************
%% Combining the model
% combine the final SIRD model into x
x = [y;y1;y2;y3];
% combien the final model based cumulative cases into case_model
cases_model = [new_cases_model,new_cases_model2,new_cases_model3,new_cases_model4];
measured = [measured_infected;measured_deaths].';
% plot the figures
figure
plot(x(:,4),'LineWidth',2);
hold on;
plot(measured_deaths,'LineWidth',2);
legend('deceased (model)','measured deceased','location','northwest')
title("comparison of cumulative decease")
figure
plot(cases_model,'LineWidth',2);
hold on;
plot(measured_infected,'LineWidth',2);
legend('cases (model)','measured infected','location','northwest')
title("comparison of cumulative infection")
% *******************************************************************
%% Constructing policies
% calculate the original cases and deaths over the omicron period
netOcase1 = new_cases_model3(21)-new_cases_model3(1);
netOdeath1 = y2(21,4)-y2(1,4);

% New update matrix from the policy
infectious_rate = 0.0032;
immune_rate = 0.25;
death_rate = 0.005;
recover_rate = 0.4;
rein_rate = 0.05;

A = [1-infectious_rate  recover_rate                            rein_rate   0;
     infectious_rate    1-(recover_rate+immune_rate+death_rate) 0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];
% run the same model but with different matrix A.
t = (85:105)-84;
x0 = A*y1(size(y1,1),:).';
y2 = model(A,B,t,x0);
new_cases_model3 = new_cases_model2(length(new_cases_model2))+cumsum((y2(:,1)*infectious_rate).');
% caculate the new cases and deaths over the omicron period.
netOcase2 = new_cases_model3(21)-new_cases_model3(1);
netOdeath2 = y2(21,4)-y2(1,4);
% calculate the percentage of change.
percent1 = (netOcase2-netOcase1)/netOcase1;
percent2 = (netOdeath2-netOdeath1)/netOdeath1;
% *******************************************************************
%% Define functions
% function that return the SIRD model using ss and lsim function.
function Y = model(A,B,t,x0)
size1 = size(t,2);
sys_sir_base = ss(A,B,eye(4),zeros(4,1),1);
Y = lsim(sys_sir_base,zeros(size1,1),t,x0);
end