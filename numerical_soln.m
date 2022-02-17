%% Defining Constants
clc; 
clear;
R_E = 6371.23 * 1000; %m
g = 9.8066; %m/s^2
rho_0 = 1.5; %kg/m^3
beta = 1/6900; %m^-1
mass = 40; %kg
% S = 2402 * 1e-6; %m^2
max_R = 765/1000; %m
S = pi * max_R * max_R; 
C_D = [0.8 1.4]; 
k_d = mass./(S * C_D); 
z_0 = 120 * 1000; %m
v_0 = 7396.6; %m/s
gamma_0 = [-1.4 -10 -40 -60]; %degrees

%% Diff Eq
V = {}; 
gamma = {}; 
Z = {}; 

tspan = [0 500];

i = 1; 

for gamma_itr = gamma_0
    y0 = [v_0 deg2rad(gamma_itr) z_0];
    [t,y] = ode45(@(t,y) odefcn(t,y,beta,k_d(1),g,R_E,rho_0), tspan, y0);

    V{i} = y(:,1); 
    gamma{i} = y(:,2); 
    Z{i} = y(:,3);

    i = i+1; 

    y0 = [v_0 deg2rad(gamma_itr) z_0];
    [t,y] = ode45(@(t,y) odefcn(t,y,beta,k_d(2),g,R_E,rho_0), tspan, y0);

    V{i} = y(:,1); 
    gamma{i} = y(:,2); 
    Z{i} = y(:,3);

    i = i+1; 

end

for i = 1:length(V)
    index = length(V{i}); 
    for j = 1:length(V{i})
        if Z{i}(j) < 0
            index = j-1; 
            break;
        end
    end
    V{i} = V{i}(1:index);
    gamma{i} = gamma{i}(1:index);
    Z{i} = Z{i}(1:index);

end
%% Plotting V(z)
figure(1)
plot(V{1}, Z{1}, '-.');
hold on; 
plot(V{2}, Z{2}, '-o');
plot(V{3}, Z{3}, '-.');
plot(V{4}, Z{4}, '-o');
plot(V{5}, Z{5}, '-.'); 
plot(V{6}, Z{6}, '-o');
plot(V{7}, Z{7}, '-.'); 
plot(V{8}, Z{8}, '-o');
legend('\gamma = -1.4 C_D = 0.8', '\gamma = -1.4 C_D = 1.4', '\gamma = -10 C_D = 0.8', '\gamma = -10 C_D = 1.4', '\gamma = -40 C_D = 0.8', '\gamma = -40 C_D = 1.4', '\gamma = -60 C_D = 0.8', '\gamma = -60 C_D = 1.4');
xlabel('Velocity (m/s)');
ylabel('Altitude (m)');
title('Flight Path');

%% Comparison between analytical and numerical solution
% case for gamma = -1.4 and C_D = 0.8

Z_analytical = linspace(120*1000, 0, 1000); 
V_analytical = v_0 * exp((0.5/(k_d(1) * beta * sin(gamma_0(1)))) * rho_0 * exp((-1 * Z_analytical * beta))); 

figure(2)
plot(V_analytical, Z_analytical, '-.'); 
hold on; 
plot(V{1}, Z{1}, '-o'); 
xlabel('Velocity (m/s)');
ylabel('Altitude (m)');
title('Flight Path Comparison for C_D = 0.8 and \gamma_0 = -1.4 degrees');
legend('Analytical', 'Numerical');

% case for gamma = -1.4 and C_D = 1.4

figure(3)
Z_analytical = linspace(120*1000, 0, 1000); 
V_analytical = v_0 * exp((0.5/(k_d(2) * beta * sin(gamma_0(1)))) * rho_0 * exp((-1 * Z_analytical * beta))); 

plot(V_analytical, Z_analytical, '-.'); 
hold on; 
plot(V{2}, Z{2}, '-o'); 
xlabel('Velocity (m/s)');
ylabel('Altitude (m)');
title('Flight Path Comparison for C_D = 1.4 and \gamma_0 = -1.4 degrees');
legend('Analytical', 'Numerical');

%% Calculating and Plotting Dynamic Pressure 
q = {};

for i = 1:length(V)
    rho_i = rho_0 * exp(-1*beta*Z{i}); 
    q{i} = 0.5 * rho_i .* (V{i}).^2;
end

figure(4)
plot(q{1}, Z{1}, '-.');
hold on; 
plot(q{2}, Z{2}, '-o');
plot(q{3}, Z{3}, '-.');
plot(q{4}, Z{4}, '-o');
plot(q{5}, Z{5}, '-.'); 
plot(q{6}, Z{6}, '-o');
plot(q{7}, Z{7}, '-.'); 
plot(q{8}, Z{8}, '-o');
legend('\gamma = -1.4 C_D = 0.8', '\gamma = -1.4 C_D = 1.4', '\gamma = -10 C_D = 0.8', '\gamma = -10 C_D = 1.4', '\gamma = -40 C_D = 0.8', '\gamma = -40 C_D = 1.4', '\gamma = -60 C_D = 0.8', '\gamma = -60 C_D = 1.4');
xlabel('Dynamic Pressure (kg/ms^2)');
ylabel('Altitude (m)');
title('Altitude v/s Dynamic Pressure');

%% Calculating Mach
[T, a, P, rho] = atmoscoesa(Z_analytical);
Mach = V_analytical./a;
figure(5)
plot(Mach, Z_analytical, '-.')
legend(' ');
xlabel('Mach Number');
ylabel('Altitude (m)');
title('Altitude v/s Mach Number');

v_table = [V_analytical(1),V_analytical(100),V_analytical(200),V_analytical(300),V_analytical(400),V_analytical(500),V_analytical(600),V_analytical(700),V_analytical(800),V_analytical(900),V_analytical(1000)].';
t_table = [T(1),T(100),T(200),T(300),T(400),T(500),T(600),T(700),T(800),T(900),T(1000)].';
rho_table = [rho(1),rho(100),rho(200),rho(300),rho(400),rho(500),rho(600),rho(700),rho(800),rho(900),rho(1000)].';
M_table = [Mach(1),Mach(100),Mach(200),Mach(300),Mach(400),Mach(500),Mach(600),Mach(700),Mach(800),Mach(900),Mach(1000)].';
z_table = [Z_analytical(1),Z_analytical(100),Z_analytical(200),Z_analytical(300),Z_analytical(400),Z_analytical(500),Z_analytical(600),Z_analytical(700),Z_analytical(800),Z_analytical(900),Z_analytical(1000)].';

table(v_table,t_table,rho_table,M_table,z_table);


%% Function
% Use V, gamma, z
function dydt = odefcn(t,y,beta,k_d,g,R_E,rho_0)
  dydt = zeros(3,1);
  dydt(1) = (-0.5*rho_0*exp(-1*beta*y(3)) * y(1)^2)/k_d - g*(1 - (y(1)^2)/(R_E*g*(1 + y(3)/R_E)))*sin(y(2));
  dydt(2) = (-1*g*(1 - (y(1)^2)/(R_E*g*(1+y(3)/R_E)))*cos(y(2)))/y(1);
  dydt(3) = y(1) * sin(y(2)); 
end 
