clc; 
clear;

% Angle of attack and entry angle
gamma_0 = -1.4;
alpha_iter = [0, 5, 10, 15, 20]; %degrees

% File to read in (capsule geometry)
file = 'CAD_capsule_3.stl';

C_D_plot = {};
C_L_plot = {};
M_plot = {};

k = 1;

for alpha = alpha_iter
    
C_D = 1.4*ones(1,500)';
C_L = zeros(1,500)';
time = linspace(0, 500, 500)';
plotting = false;

disp('Solving Diff Eq'); 
[M1, V, time, rho, P] = convergence_numsoln(C_D, C_L, time, gamma_0);

% clipping Mach number
clip_index = length(M1); 
for i = 1:length(M1)
    if M1(i) < 2
        clip_index = i-1; 
        break
    end
end

disp('Calculating Coefficients')
[C_D, C_L, pressures, areas, f_z, f_x, Cps] = pressure_calc(M1, V, alpha, plotting, rho, P, file);

C_D = C_D(1:clip_index); 
C_L = C_L(1:clip_index);
M1 = M1(1:clip_index);

C_D_plot{k} = C_D;
C_L_plot{k} = C_L;
M_plot{k} = M1;
disp(k)
k = k+1;
end

