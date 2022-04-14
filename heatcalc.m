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
%% Calculate Dynamic Pressure
Dynamic_Pressure = 0.5 * rho_traj.* V.^2;
%% Calculate Load Factor
acc = zeros((length(V)-1),1);
Z_acc = zeros((length(V)-1),1);
for j = 1:(length(V)-1)
   acc(j) = (V(j)-V(j+1))/(time(j)-time(j+1));
   Z_acc(j) = (Z(j)+Z(j+1))/2;
end
g = -9.81;
Load_Factor = acc/g;
%% Plot optimized V vs Z
figure(1)
set(gcf,'color','w');
plot(Mach, Z)
hold on
plot(Load_Factor, Z_acc)
title('Trajectory Properties v/s Altitude (Optimized Flight Path')
legend('Mach', 'Load Factor')
xlabel('Property')
ylabel('Altitude (m)')

%% Find Cp_0 along the Flight Path
gamma = gamma_atmos;
M1 = Mach;
% nose pressure, contingent on shock theory and M1 dependance
P2_P1 = 1 + 2*gamma/(gamma+1)*(M1.^2-1);
Pstag_P2 = ((gamma+1)^2*M1.^2./(4*gamma*M1.^2-2*(gamma-1))).^3;
C_p0 = 2./(gamma*M1.^2).*(P2_P1 .* Pstag_P2-1);

%% Calculate Stagnation Pressure along Flight Path

R_n = 0.272; %m
% c_p = 1.00; % kJ/kg
N = 0.5;
M = 3;
r = 0.71^0.5; %laminar flow where 0.71 is Prandtl as r = Pr^n
T_e = [800 900 1000 1100 1200]; %K
q_w = {};

for i = 1:5
    % Fix wall temperature
    T_w = T_e(i);
    % Tauber Menees
    %T_aw = T_traj + r*V.^2./(2*C_p0);
    T_aw = 1800;
    C = (1.83*10^-8)*R_n^-0.5*(1-T_w/T_aw)*10000;
    q_w{i} = (C*rho_traj.^N.*V.^M)/1000;
end

%% Calculate Blackbody Radiation along Traj

eps = 0.7;
sigma = 5.6695 * 10^-8; % Stefan Boltzmann [W/m^2/K^4]
q_r = {};
for i = 1:5
    % Fix wall temperature
    T_w = T_e(i);
    % Tauber Menees
    q_r{i} = (eps*sigma*T_w^4 * ones(length(Z),1))/1000;
end

%% Correct q_total

q_total = {};
for i = 1:5
    q_total{i} = (q_w{i} - q_r{i});
end

%% Plotting q_w and q_total

figure(2)
set(gcf,'color','w');
plot(q_w{1}, Z)
hold on
plot(Dynamic_Pressure, Z)
title('Trajectory Properties v/s Altitude')
xlabel('Properties')
ylabel('Altitude (m)')
legend('Q_{wall}','Dynamic Pressure')

%% Plotting q_total across T

figure(4)
plot(q_total{1}, Z)
hold on
plot(q_total{2}, Z)
hold on
plot(q_total{3}, Z)
hold on
plot(q_total{4}, Z)
hold on
plot(q_total{5}, Z)
hold off

title('Stagnation Flux at Various Temperatures along Trajectory')
xlabel('q_w (W/cm^2)')
ylabel('Altitude (m)')
legend('800','900','1000','1100','1200')

%% Sutton and Graves Stagnation Heat Exchange
q_w_c = (1.1415*10^-4).*(rho_traj.^(1/2)).*V.^3./R_n^(1/2)/10000;
epsilon =  0.9;
K1 = 372.6;
K2 = 8.5;
K3 = 1.6;
q_w_r_stag = (R_n*K1*((3.28084*10^-4).*V).^K2).*((rho_traj./1.22499915588771).^K3)*10000;
T_w_sutgrav = (((q_w_c + q_w_r_stag)/(epsilon*sigma)).^0.25);
q_r_wall = (epsilon*sigma*T_w_sutgrav.^4);
%q_w_cond = ((q_w_c + q_w_r_stag - q_r_wall));
q_w_cond = 0; %assumption

figure(5)
plot(q_w_c, Z)
hold on
plot(q_w_r_stag, Z)
plot(q_r_wall, Z)
%plot(q_w_cond, Z)
hold off
title('Heat Exchange at Stagnation Point along the Trajectory')
xlabel('qs (kW/m^2)')
ylabel('Altitude (m)')
legend('qc', 'qr', 'qwall', 'qcond')

%% Tauber-Menees Heat Rate Along Body
N = 0.5;
M = 3.2;
% used CAD for measurements
r_cone_1 = 0.270356; %m
r_cone_2 = 0.0692;
S_nose = 0.3437976; %estimate
S_cone = 0.698695; %m rest of S after 
r = 0.71^0.5; %laminar flow
A = 0.664; % (8.53.1)
m = 0.5;

%solve for Taw and Tw
%T_aw = T_traj.*(1 + r.*((gamma_atmos - 1)/2).*Mach.^2);
T_aw = T_traj + r*V.^2./(2*C_p0);
%T_aw = 100000;
%T2_Ttraj = ((2*gamma_atmos * (gamma_atmos - 1))/(gamma_atmos + 1)^2).*Mach;
%Ttraj_T2 = T2_Ttraj.^-1;
% T2 = T2_Ttraj.*T_traj;
% M2 = ((((gamma_atmos - 1)*Mach.^2 + 2))./(2*gamma_atmos.*Mach.^2-(gamma_atmos - 1))).^0.5;
%T_w = T2.*((Ttraj_T2 - 0.16*r*((gamma_atmos - 1)/2).*M2.^2 - 1))/0.55;
%T_w = [100 500 1000 1500 2000]; %K
T_w = 1000;
theta_wedge = pi/4;

% discretizing S along the cone
S_cone_inc = [];
for i = 1:200
    S_cone_inc(i) = S_nose + i*S_cone/200;
end

C_cyl = [];
for i = 1:length(S_cone_inc)
    for j = 1:length(rho_traj)
        C_cyl(j,i) = (2.42*10^-9)*(cos(theta_wedge)^0.5)*sin(theta_wedge)*(S_cone_inc(i)^(-0.5)).*(1-T_w./T_aw(j));
    end
end

q_w_cone = [];
for i = 1:length(S_cone_inc)
    for j = 1:length(rho_traj)
        q_w_cone(j,i) =  (C_cyl(j,i)*10)*rho_traj(j)^N*V(j)^M /10000;
        
    end
end
q_w_cones = sum(q_w_cone, 2); % getting total q on cone

figure(6)
plot(q_w_cones(1:12000), Z(1:12000))
title('Heat Exchange on Cone along the Trajectory')
xlabel('qwcone (W/cm^2)')
ylabel('Altitude (m)')
