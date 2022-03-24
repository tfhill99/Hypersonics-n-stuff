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

%% Calculate Mach, Rho, P_inf along traj
% First need T from atmo function
gamma_atmos = 1.4;
R_atmos = 287; %J/kg/K

[Z_total, Z_L, Z_U, T, P, rho, c, g_atmos, mu, nu, k, n, n_sum]  = atmo(120, 0.01, 1);

% flip trajectory
Z_total = flip(Z_total)*1000; %convert to m
P = flip(P);
rho = flip(rho);

T_traj = interp1(Z_total, T, Z);
P_traj = interp1(Z_total, P, Z);
rho_traj = interp1(Z_total, rho, Z);

c_traj = sqrt(gamma_atmos*R_atmos*T_traj);

Mach = V./c_traj;
%% Plot optimized V vs Z

figure(1)
plot(Z, V)
title('Velocity v/s Altitude (Optimized Flight Path)')
xlabel('Altitude (m)')
ylabel('Velocity (m/s)')
legend('Velocity')

%% Find Cp along the Flight Path

alpha_init = rad2deg(alpha_r(1));
plotting = false;
file = 'CAD_capsule_3.stl';

[C_D, C_L, C_M, pressures, areas, f_z, f_x, Cps, C_p0] = pressure_calc(Mach, V, alpha_init, plotting, rho, P, file);

%% Calculate Stagnation Pressure along Flight Path

R_n = 0.272; %m
c_p = 1.00; % kJ/kg
N = 0.5;
M = 3;
r = 0.71^0.5; %laminar flow where 0.71 is Prandtl as r = Pr^n
T_e = [100 500 1000 1500 2000]; %K
q_w = {};

for i = 1:5
    % Fix wall temperature
    T_w = T_e(i);
    % Tauber Menees
    T_aw = T_w + r*V.^2/(2*c_p);
    C = 1.29*10^-4*R_n^-0.5*(1-T_w/T_aw);
    q_w{i} = C*rho_traj.^N.*V.^M;
end

%% Calculate Blackbody Radiation along Traj

eps = 0.7;
sigma = 5.6695 * 10^-8; % Stefan Boltzmann [W/m^2/K^4]
q_r = {};
for i = 1:5
    % Fix wall temperature
    T_w = T_e(i);
    % Tauber Menees
    q_r{i} = eps*sigma*T_w^4 * ones(length(Z),1);
end

%% Correct q_total

q_total = {};
for i = 1:5
    q_total{i} = q_w{i} - q_r{i};
end

%% Plotting q_w and q_total

figure(2)
plot(Z, q_w{1})
hold on
plot(Z, q_total{1})
title('Stagnation Flux along Trajectory, T_w = 800 K')
xlabel('Altitude (m)')
ylabel('q_w (W/cm^2)')
legend('Without Blackbody Radiation','With Blackbody Radiation')

%% Plotting q_total across T

figure(3)
plot(Z, q_total{1})
hold on
plot(Z, q_total{2})
hold on
plot(Z, q_total{3})
hold on
plot(Z, q_total{4})
hold on
plot(Z, q_total{5})
hold on
title('Stagnation Flux at Various Temperatures along Trajectory')
xlabel('Altitude (m)')
ylabel('q_w (W/cm^2)')
legend('800','900','1000','1100','1200')