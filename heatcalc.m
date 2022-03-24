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

R_n = 2.72; %m
c_p = 1.00; % kJ/kg
T_e = [800 900 1000 1100 1200]';
q_w = {};

for i = 1:5
    T_w = T_e(i);
    % Tauber Menees
    N = 0.5;
    M = 3;
    r = 0.71^0.5; %laminar flow where 0.71 is Prandtl as r = Pr^n
    T_aw = T_w + r*V.^2/(2*c_p);
    C = 1.29*10^-4*R_n^-0.5*(1-T_w/T_aw);
    q_w{i} = C*rho_traj.^N.*V.^M;
end

%% Plotting q

figure(2)
plot(Z, q_w{1})
hold on
plot(Z, q_w{2})
hold on
plot(Z, q_w{3})
hold on
plot(Z, q_w{4})
hold on
plot(Z, q_w{5})
hold on
title('Stagnation Heating at Various Temperatures along Trajectory')
xlabel('Altitude (m)')
ylabel('q_w')
legend('800','900','1000','1100','1200')