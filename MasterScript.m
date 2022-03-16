clc; 
clear;
% Initial ballistic C_D and C_L
C_D = 1.4*ones(1,500)';
C_L = zeros(1,500)';
% Angle of attack and entry angle
gamma_0 = -1.4;
alpha = 10;
time = linspace(0, 500, 500)';
% Arrays for tracking C_D and C_L convergence
C_D_change = {};
C_D_change{1} = C_D;
C_L_change = {};
C_L_change{1} = C_L;
plotting = false;
% File to read in (capsule geometry)
file = 'CAD_capsule_3.stl';
[M1, V, time, rho, P] = convergence_numsoln(C_D, C_L, time, gamma_0);
disp('Step 1 done yay')
[C_D, C_L] = pressure_calc(M1, V, alpha, plotting, rho, P, file);
disp('BC')
C_D_change{2} = C_D;
C_L_change{2} = C_L;

%{
i = 2;
while C_D_change{i} - C_D_change{i-1} > 1
    [M1, V, time, rho, P] = convergence_numsoln(C_D, C_L, time, gamma_0);
    [C_D, C_L] = pressure_calc(M1, V, alpha, plotting, rho, P, file);
    C_D_change{i+1} = C_D;
    C_L_change{i+1} = C_L;
    i = i + 1;
    i
end
%}


