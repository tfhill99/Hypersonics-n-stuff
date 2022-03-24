clear all
clc

%% Read in variables
traj = table2array(readtable('OptimizedPathDatabase.xlsx'));
time = traj(3:end,1);
V = traj(3:end,2);
alpha_r = traj(3:end,3);
q = traj(3:end,4);
X = traj(3:end,5);
Z = traj(3:end,6);
theta_r = traj(3:end,7);
gamma_r = traj(3:end,8);

%% Plot optimized V vs Z

figure(1)
plot(Z, V)
title('Velocity v/s Altitude (Optimized Flight Path)')
xlabel('Altitude (m)')
ylabel('Velocity (m/s)')
legend('Velocity')

%% Calculate Stagnation Pressure along Flight Path

% Tauber Menees
N = 0.5;
M = 3;

q_w = C*rho_inf.^N*V_inf.^M;