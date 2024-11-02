clearvars;
close all;

% Question 1: Simulate the model described in Section 9.3 of the book
X = zeros(4,100);
X(:,1) = [1,0,0,0]';
t = 1:200;
A = [0.95,0.04,0,0; 
    0.05,0.85,0,0;
    0,0.10,1,0;
    0,0.01,0,1];

for i = 2:200
    X(:,i) = A * X(:,i-1);
end

hold on;
plot(t, X(1,:),'LineWidth',2);
plot(t, X(2,:),'LineWidth',2);
plot(t, X(3,:),'LineWidth',2);
plot(t, X(4,:),'LineWidth',2);
xlabel("time");
ylabel("x_t");
legend(["Susceptible", "Infection", "Recovered", "Deceased"]);
title("Base SIRD Model")
fprintf("After about 100 days, the state converges to 1 and 10% of the population deceased.");

% Question 2: Implelement the senario in which 1 percent of the recovered
% population becomes susceptible 
A = [0.95, 0.04,  0.01, 0; 
     0.05, 0.85,  0,   0;
     0,    0.10,  0.99, 0;
     0,    0.01,  0,   1];
t = 1:1000;

for i = 2:1000
    X(:,i) = A * X(:,i-1);
end


figure;
hold on;
plot(t, X(1,:),'LineWidth',2);
plot(t, X(2,:),'LineWidth',2);
plot(t, X(3,:),'LineWidth',2);
plot(t, X(4,:),'LineWidth',2);
xlabel("time");
ylabel("x_t");
legend(["Susceptible", "Infection", "Recovered", "Deceased"]);
title("SIRD mode with re-infection")

% Question 3
% Below is the code in the base_sir.m file, which uses ss and lsim function
% to simulate
A = [0.95 0.04 0 0; 0.05 0.85 0 0; 0 0.1 1 0; 0 0.01 0 1];

% The following matrix is needed to use the lsim function to simulate the
% system in question
B = zeros(4,1);

% initial conditions (i.e., values of S, I, R, D at t=0).
x0 = [0.9 0.1 0 0];

% Here is a compact way to simulate a linear dynamical system.
% Type 'help ss', 'help lsim', etc., to learn about how these functions work!!
sys_sir_base = ss(A,B,eye(4),zeros(4,1),1);
Y = lsim(sys_sir_base,zeros(200,1),linspace(0,199,200),x0);

% plot the output trajectory
figure;
plot(Y);
legend('S','I','R','D');
xlabel('Time')
ylabel('Percentage Population');

%The output using the ss and lsim functions has the same outcome as our own in part i.