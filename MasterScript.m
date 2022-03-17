clc; 
clear;

% Angle of attack and entry angle
gamma_0 = -1.4;
alpha = 0; %degrees

% Initial ballistic C_D and C_L
C_D = 1.4*ones(1,500)';
C_L = zeros(1,500)';
time = linspace(0, 500, 500)';

% Arrays for tracking C_D and C_L convergence
C_D_change = {};
C_D_change{1} = C_D;
C_L_change = {};
C_L_change{1} = C_L;
plotting = false;

% File to read in (capsule geometry)
file = 'CAD_capsule_3.stl';

disp('Solving Diff Eq'); 
[M1, V, Z, time, rho, P, T, A] = convergence_numsoln(C_D, C_L, time, gamma_0);

% clipping Mach number
clip_index = length(M1); 
for i = 1:length(M1)
    if M1(i) < 2
        clip_index = i-1; 
        break
    end
end

%{
M1 = M1(1:clip_index); 
V = V(1:clip_index); 
time = time(1:clip_index); 
rho = rho(1:clip_index); 
P = P(1:clip_index); 
%}

disp('Calculating Coefficients')
[C_D, C_L, pressures, areas, f_z, f_x, Cps] = pressure_calc(M1, V, alpha, plotting, rho, P, file);
C_D_change{2} = C_D;
C_L_change{2} = C_L;

i = 2; 

A = C_D_change{i}; 
B = C_D_change{i-1}; 
A = A(1:clip_index); 
B = B(1:clip_index); 

while any(A - B > 0.05)
    disp('Solving Diff Eq'); 
    [M1, V, Z, time, rho, P, T, A] = convergence_numsoln(C_D, C_L, time, gamma_0);

    % clipping Mach number
    
    clip_index = length(M1); 
    for j = 1:length(M1)
        if M1(j) < 2
            clip_index = j-1; 
            break
        end
    end

    disp('Calculating Coefficients')
    [C_D, C_L, pressures, areas, f_z, f_x, Cps] = pressure_calc(M1, V, alpha, plotting, rho, P, file);
    C_D_change{i+1} = C_D;
    C_L_change{i+1} = C_L;
    i = i + 1;
    disp(i)
    A = C_D_change{i}; 
    B = C_D_change{i-1};
    A = A(1:clip_index); 
    B = B(1:clip_index); 
end

C_D_final = C_D_change{i}; 
C_D_final = C_D_final(1:clip_index); 

C_L_final = C_L_change{i}; 
C_L_final = C_L_final(1:clip_index);

M1 = M1(1:clip_index);
disp('done')

%% Plotting Final Trajectory

figure(1)
set(gcf,'color','w');
plot(M1, C_D_final);
xlabel('Mach Number');
ylabel('C_d');
title('Converged C_D versus Mach for AoA = 0');

%% Test Matrix

% Altitude, Velocity, Mach, Acceleration, Temperature, Density



