clc;
clearvars;
close all;
% *******************************************************************
load("COVID_STL.mat");
POP_2 = 1901723;
days = 300;
t = 1:days;

% ********************************************************************
%% 4.2.1 Part III, Question 1
X = zeros(4,days);
X(:,1) = [POP_STL,0,0,0]';
Xt = zeros(4,days);
Xt(:,1) = [POP_STL,0,0,0]';
X2 = zeros(4,days);
X2(:,1) = [POP_2,0,0,0]';
Xt2 = zeros(4,days);
Xt2(:,1) = [POP_2,0,0,0]';

% Set up A1 with arbitrary rates
infectious_rate = 0.003;
immune_rate = 0.092;
death_rate = 0.0086;
recover_rate = 0.20;
rein_rate = 0.005;
A1 = [1-infectious_rate recover_rate                            rein_rate   0;
     infectious_rate    1-(recover_rate+immune_rate+death_rate) 0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];

% Set up A2 with arbitrary rates
infectious_rate_2 = 0.002;
immune_rate_2 = 0.0921;
death_rate_2 = 0.018;
recover_rate_2 = 0.201;
rein_rate_2 = 0.003;
A2 = [1-infectious_rate_2  recover_rate_2                                rein_rate_2   0;
     infectious_rate_2     1-(recover_rate_2+immune_rate_2+death_rate_2) 0             0;
     0                     immune_rate_2                                 1-rein_rate_2 0;
     0                     death_rate_2                                  0             1];

% Randomly generate the number of people leaving each city
u1 = ones(size(t,2),1)*POP_STL;
u2 = ones(size(t,2),1)*POP_2;
for i = 1:size(u1,1)
    u1(i) = u1(i) * (0.001 + (0.008-0.001) * rand);
    u2(i) = u2(i) * (0.001 + (0.008-0.001) * rand);
end

% Update population for each time interval
newPop_stl = zeros(size(t));
newPop_2 = newPop_stl;
newPop_stl(:,1) = POP_STL;
newPop_2(:,1) = POP_2;

for i = 1:days-1
    %% Without travel SIRD number
    X(:,i+1) = A1 * X(:,i);
    X2(:,i+1) = A2 * X2(:,i);
    
    %% With travel
    % new total population
    newPop_stl(:,i+1) = newPop_stl(:,i) - u1(i,:) + u2(i,:);
    newPop_2(:,i+1) = newPop_2(:,i) - u2(i,:) + u1(i,:);

    % New SIRD number
    B1 = [Xt(1,i)/newPop_stl(i); Xt(2,i)/newPop_stl(i); Xt(3,i)/newPop_stl(i); 0];
    B2 = [Xt2(1,i)/newPop_2(i); Xt2(2,i)/newPop_2(i); Xt2(3,i)/newPop_2(i); 0];
    leaving_stl = B1 * u1(i,:);
    leaving_2 = B2 * u2(i,:);
    Xt(:,i+1) = A1 * Xt(:,i) - leaving_stl + leaving_2;
    Xt2(:,i+1) = A2 * Xt2(:,i) - leaving_2 + leaving_stl;
end

% Calculate the cumulative cases number proportional to the total population 
Cumulative_Case = cumsum(X(1,:)*infectious_rate');
Cumulative_Case2 = cumsum(X2(1,:)*infectious_rate_2');
Cumulative_Case_with_Travel = cumsum(Xt(1,:)*infectious_rate');
Cumulative_Case_with_Travel2 = cumsum(Xt2(1,:)*infectious_rate_2');

% Calculate the cumulative cases for the combined of the two cities.
Total = (Cumulative_Case + Cumulative_Case2)/(POP_STL+POP_2);
Total_Travel = (Cumulative_Case_with_Travel+Cumulative_Case_with_Travel2)/(POP_STL+POP_2);

Cumulative_Case = Cumulative_Case/POP_STL;
Cumulative_Case2 = Cumulative_Case2/POP_2;

% Calculate the cumulative cases number proportional to the population each
% day
for i = 1:days
    Cumulative_Case_with_Travel(:,i) = Cumulative_Case_with_Travel(:,i)/newPop_stl(:,i);
    Cumulative_Case_with_Travel2(:,i) = Cumulative_Case_with_Travel2(:,i)/newPop_2(:,i);
end

% Plot the data
hold on;
plot(t,Cumulative_Case);
plot(t,Cumulative_Case_with_Travel);
plot(t,Cumulative_Case2);
plot(t,Cumulative_Case_with_Travel2);
title("Cumulative Case Proportion of Interaction of two Distinct Population")
xlabel("time");
ylabel("Cumulative Case proportional to the population");
legend("infected",'infected with travel',"infected 2",'infected with travel 2');

figure;
hold on;
plot(t,Total);
plot(t,Total_Travel);
title("Total Cumulative Case Proportions")
xlabel("time");
ylabel("Cumulative Case proportional to the population");
legend("total cumulative", 'total cumulative with travel');

% ********************************************************************
%% 4.2.1 Part III, Question 2
% setup the time interval for two phases
section1 = 150;
seciton2 = days-section1;
t1 = 1 : section1;
t2 = section1+1 : days;

% generate the number of people leaving the city after a more strict travel
% policy
for i = section1+1 : size(u1,1)
    u1(i) = u1(i) * (0.0005 + (0.005-0.0005) * rand);
    u2(i) = u2(i) * (0.0005 + (0.005-0.0005) * rand);
end

for i = 1:days-1
    %% With travel
    % new total population
    newPop_stl(:,i+1) = newPop_stl(:,i) - u1(i,:) + u2(i,:);
    newPop_2(:,i+1) = newPop_2(:,i) - u2(i,:) + u1(i,:);

    % New SIRD number
    B1 = [Xt(1,i)/newPop_stl(i); 0; Xt(3,i)/newPop_stl(i); 0];
    B2 = [Xt2(1,i)/newPop_2(i); 0; Xt2(3,i)/newPop_2(i); 0];
    leaving_stl = B1 * u1(i,:);
    leaving_2 = B2 * u2(i,:);
    Xt(:,i+1) = A1 * Xt(:,i) - leaving_stl + leaving_2;
    Xt2(:,i+1) = A2 * Xt2(:,i) - leaving_2 + leaving_stl;
end

% Calculate the cumulative cases number proportional to the total population 
Cumulative_Case_with_Travel = cumsum(Xt(1,:)*infectious_rate');
Cumulative_Case_with_Travel2 = cumsum(Xt2(1,:)*infectious_rate_2');

% Calculate the cumulative cases for the combined of the two cities.
Total_Travel = (Cumulative_Case_with_Travel+Cumulative_Case_with_Travel2)/(POP_STL+POP_2);

% Calculate the cumulative cases number proportional to the population each
% day
for i = 1:days
    Cumulative_Case_with_Travel(:,i) = Cumulative_Case_with_Travel(:,i)/newPop_stl(:,i);
    Cumulative_Case_with_Travel2(:,i) = Cumulative_Case_with_Travel2(:,i)/newPop_2(:,i);
end


% Plot the data
figure;
hold on;
plot(t,Cumulative_Case);
plot(t,Cumulative_Case_with_Travel);
plot(t,Cumulative_Case2);
plot(t,Cumulative_Case_with_Travel2);
title("Cumulative Case Proportion of Interaction of two Distinct Population")
xlabel("time");
ylabel("Cumulative Case proportional to the population");
legend("infected",'infected with travel',"infected 2",'infected with travel 2');
snapnow

figure;
hold on;
plot(t,Total);
plot(t,Total_Travel);
title("Total Cumulative Case Proportions")
xlabel("time");
ylabel("Cumulative Case proportional to the population");
legend("total cumulative", 'total cumulative with travel');
snapnow
% *********************************************************************
