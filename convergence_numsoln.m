function [Mach, V, Z, time, rho, P, T, A, Re] = convergence_numsoln(C_D, C_L, time, gamma_0)
%% Defining Constants
R_E = 6371.23 * 1000; %m
g = 9.8066; %m/s^2
rho_0 = 1.5; %kg/m^3
beta = 1/6900; %m^-1
mass = 40; %kg
max_R = 765/1000; %m
R = 287; % J/kg*K
gamma_atmos = 1.4; % constant for air
S = pi * max_R * max_R; 
z_0 = 120 * 1000; %m
v_0 = 7396.6; %m/s
kdt = time;
klt = time;

%% Diff Eq
k_d = (mass/S)./C_D;
k_l = (mass/S)./C_L;

tspan = [0 500];
t0 = linspace(0, 300, 300);

y0 = [v_0 deg2rad(gamma_0) z_0];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t,y] = ode45(@(t,y) myode(t,y,k_d,kdt,k_l,klt,beta,g,R_E,rho_0), tspan, y0);
     
V = interp1(t, y(:,1),t0)';
gamma = interp1(t, y(:,2),t0)';
Z = interp1(t, y(:,3),t0)';
time = t0';


%% Mach
[Z_total, Z_L, Z_U, T, P, rho, c, g_atmos, mu, nu, k, n, n_sum]  = atmo(120, 0.01, 1);
% Outputs go from index(1) = sea level to index(12001) = 120,000m

Z_L = flip(Z_L)*1000;
c = flip(c);
index = 0;
% Find which of the numerical Z values are in range of atmos 86,000m 
for i = 1:length(Z)
   if Z(i) > Z_L(1)
       index = index + 1;
   end
end
Z_inlimit = Z(index+1:end);
Z_extend = Z(1:index);

% Interpolate the respective speed of sound values for numerical solution
c_inlimit_interp = interp1(Z_L, c, Z_inlimit);


% Find the temperatures for 86000 - 120000 km
T_extend = T(length(Z_L):end);
T_extend = flip(T_extend);

% Extrapolate the speed of sound values from 86000-120000
c_extend = sqrt(gamma_atmos*R*T_extend);

% Interpolate the respective speed of sounds values for numerical solution
Z_U = flip(Z_U)*1000;
c_extend_interp = interp1(Z_U, c_extend, Z_extend);

% Combine the two for the speed of sound along the whole trajectory 
c_combined = cat(1,c_extend_interp,c_inlimit_interp);

% Calculate Mach number for the whole trajectory
Mach = (V./c_combined);

%% Interpolate rho and P and T. Needed for pressure_calc function input and final test matrix

% Ensure all variables are going from h --> 0
rho = flip(rho);
P = flip(P);
T = flip(T);
Z_total = flip(Z_total)*1000; % Scale to m

P = interp1(Z_total,P,Z);
rho = interp1(Z_total,rho,Z);
T = interp1(Z_total, T, Z);

%% Calculate Acceleration. Needed for final test matrix

A = {};
for j = 1:(length(V)-1)
   A{j} = (V(j)-V(j+1))/(time(j)-time(j+1));
end
A = A.';

%% Calculate Reynolds Number
h = 562.32; %m
V_clip = V(128:end);
rho_clip = rho(128:end);
Z_clip = Z(128:end);
mu_interp = interp1(Z_L, mu, Z_clip);

Re = h*V_clip.*rho_clip./mu_interp;
end
%% Function
function dydt = myode(t,y,k_d,kdt,k_l,klt,beta,g,R_E,rho_0)
    k_d = interp1(kdt,k_d,t); % Interpolate the data set (kd) at time t. 
    k_l = interp1(klt,k_l,t); % Interpolate the data set (kl) at time t. 
    dydt = zeros(3,1);
    dydt(1) = (-0.5*rho_0*exp(-1*beta*y(3)) * y(1)^2)/k_d - g*(1 - (y(1)^2)/(R_E*g*(1 + y(3)/R_E)))*sin(y(2));
    dydt(2) = (0.5*rho_0*exp(-1*beta*y(3)) * y(1))/k_l + (-1*g*(1 - (y(1)^2)/(R_E*g*(1+y(3)/R_E)))*cos(y(2)))/y(1);
    dydt(3) = y(1) * sin(y(2)); 
end
