clc;
clear;

% Angle of attack and entry angle
gamma_0 = -1.4;
alpha_iter = [0 10 15 20 25]; %degrees

% Arrays for tracking C_D and C_L convergence
C_D_change = {};
C_L_change = {};
C_M_change = {};
plotting = false;

% plotting arrays
C_D_plot = {};
C_L_plot = {};
C_M_plot = {};
M_plot = {};

% File to read in (capsule geometry)
file = 'CAD_capsule_3.stl';

count = 1;

% pressures_save = {}; 
% Cps_save = {}; 

for u = 1:length(alpha_iter)
C_D = 1.4*ones(1,300)';
C_L = zeros(1,300)';
C_M = zeros(1,300)';
C_D_change{1} = C_D;
C_L_change{1} = C_L;
C_M_change{1} = C_M;

time = linspace(0, 300, 300)';
alpha = alpha_iter(u);
fprintf('alpha =  %d', alpha);
disp('Solving Diff Eq');
[M1, V, Z, time, rho, P, T, Accel, Re] = convergence_numsoln(C_D, C_L, time, gamma_0);

clip_index = length(M1); 
for i = 1:length(M1)
    if M1(i) < 2
        clip_index = i-1; 
        break
    end
end

M1 = M1(1:clip_index); 
V = V(1:clip_index);  
rho = rho(1:clip_index); 
P = P(1:clip_index); 
time_clip = time(1:clip_index);
time_expand = linspace(0, time(clip_index), 300)';

M1 = interp1(time_clip, M1, time_expand);
V = interp1(time_clip, V, time_expand);
rho = interp1(time_clip, rho, time_expand);
P = interp1(time_clip, P, time_expand);

disp('Pressure Calc');
[C_D, C_L, C_M, pressures, areas, f_z, f_x, Cps] = pressure_calc(M1, V, alpha, plotting, rho, P, file);

C_D_change{2} = C_D;
C_L_change{2} = C_L;
C_M_change{2} = C_M;

k = 2; 

A = C_D_change{k};
B = C_D_change{k-1};

while any(A - B > 0.05)
    disp('Solving Diff Eq'); 
    [M1, V, Z, time, rho, P, T, Accel, Re] = convergence_numsoln(C_D, C_L, time, gamma_0);
    
    clip_index = length(M1); 
    for i = 1:length(M1)
        if M1(i) < 2
            clip_index = i-1; 
            break
        end
    end

    M1 = M1(1:clip_index); 
    V = V(1:clip_index);  
    rho = rho(1:clip_index); 
    P = P(1:clip_index); 
    time_clip = time(1:clip_index);
    time_expand = linspace(0, time(clip_index), 300)';

    M1 = interp1(time_clip, M1, time_expand);
    V = interp1(time_clip, V, time_expand);
    rho = interp1(time_clip, rho, time_expand);
    P = interp1(time_clip, P, time_expand);

    disp('Pressure Calc');
    [C_D, C_L, C_M, pressures, areas, f_z, f_x, Cps] = pressure_calc(M1, V, alpha, plotting, rho, P, file);
    %pressures_save{u} = pressures; 
    %Cps_save{u} = Cps; ]
    k = k + 1;
    C_D_change{k} = C_D;
    C_L_change{k} = C_L;
    C_M_change{k} = C_M;
    fprintf('k = %d',k);
    A = C_D_change{k};
    B = C_D_change{k-1};
    
end

C_D_final = C_D_change{k}; 
C_L_final = C_L_change{k}; 
C_M_final = C_M_change{k};

disp('Done');

C_D_plot{count} = C_D_final;
C_L_plot{count} = C_L_final;
C_M_plot{count} = C_M_final;
M_plot{count} = M1;

count = count+1;
end

%% Plotting Final Trajectory

figure(1)
set(gcf,'color','w');
plot(M1, C_D_final);
xlabel('Mach Number');
ylabel('C_d');
title('Converged C_D versus Mach for AoA = 25');

%% Test Matrix

% Altitude, Velocity, Mach, Acceleration, Temperature, Density
total = length(M1);
increment = cast(total/7, "uint8");

Z_table = [Z(1), Z(1+increment), Z(1+2*increment), Z(1+3*increment), Z(1+4*increment), Z(1+5*increment), Z(1+6*increment), Z(total)].';
v_table = [V(1), V(1+increment), V(1+2*increment), V(1+3*increment), V(1+4*increment), V(1+5*increment), V(1+6*increment), V(total)].';
M_table = [M1(1), M1(1+increment), M1(1+2*increment), M1(1+3*increment), M1(1+4*increment), M1(1+5*increment), M1(1+6*increment), M1(total)].';
acc_table = [Accel(1), Accel(1+increment), Accel(1+2*increment), Accel(1+3*increment), Accel(1+4*increment), Accel(1+5*increment), Accel(1+6*increment), Accel(total-1)].';
t_table = [T(1), T(1+increment), T(1+2*increment), T(1+3*increment), T(1+4*increment), T(1+5*increment), T(1+6*increment), T(total)].';
rho_table = [rho(1), rho(1+increment), rho(1+2*increment), rho(1+3*increment), rho(1+4*increment), rho(1+5*increment), rho(1+6*increment), rho(total)].';
C_d_table = [C_D_final(1), C_D_final(1+increment), C_D_final(1+2*increment), C_D_final(1+3*increment), C_D_final(1+4*increment), C_D_final(1+5*increment), C_D_final(1+6*increment), C_D_final(total)].';
C_l_table = [C_L_final(1), C_L_final(1+increment), C_L_final(1+2*increment), C_L_final(1+3*increment), C_L_final(1+4*increment), C_L_final(1+5*increment), C_L_final(1+6*increment), C_L_final(total)].';
C_m_table = [C_M_final(1), C_M_final(1+increment), C_M_final(1+2*increment), C_M_final(1+3*increment), C_M_final(1+4*increment), C_M_final(1+5*increment), C_M_final(1+6*increment), C_M_final(total)].';
table(Z_table,v_table,M_table,acc_table,t_table,rho_table, C_d_table, C_l_table, C_m_table)

%% Plotting
figure(1)
set(gcf,'color','w');
for i = 1:length(alpha_iter)
    plot(M_plot{i}, C_D_plot{i});
    hold on; 
end
legend('\alpha = 0', '\alpha = 5', '\alpha = 10', '\alpha = 15', '\alpha = 20');
xlabel('Mach Number');
ylabel('Drag Coefficient');
title('Flight Path');

figure(2)
set(gcf,'color','w');
for i = 1:length(alpha_iter)
    plot(M_plot{i}, C_L_plot{i});
    hold on; 
end
legend('\alpha = 0', '\alpha = 5', '\alpha = 10', '\alpha = 15', '\alpha = 20');
xlabel('Mach Number');
ylabel('Lift Coefficient');
title('Flight Path');

figure(3)
set(gcf,'color','w');
for i = 1:length(alpha_iter)
    plot(M_plot{i}, C_M_plot{i});
    hold on; 
end
legend('\alpha = 0', '\alpha = 5', '\alpha = 10', '\alpha = 15', '\alpha = 20');
xlabel('Mach Number');
ylabel('Pitching Moment Coefficient');
title('Flight Path');

%% Writing CSV Files
for i = 1:length(alpha_iter)
    D = [M_plot{i} C_D_plot{i} C_L_plot{i} C_M_plot{i}];
    name = "Mach_CD_CL__CM_at_" + alpha_iter(i); 
    csvwrite(name, D); 

    HAHA = "CpsByMachOnMesh_" + alpha_iter(i); 
    HAHA2 = "PressuresByMachOnMesh_" + alpha_iter(i);

    csvwrite(HAHA, Cps_save{i}); 
    csvwrite(HAHA2, pressures_save{i}); 
end
