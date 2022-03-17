clc; 
clear;

% Angle of attack and entry angle
gamma_iter = [-1.4, -10, -40, -60];
C_D_iter = [0.8, 1.4];
alpha = 0; %degrees

% File to read in (capsule geometry)
file = 'CAD_capsule_3.stl';

C_D_plot = {};
C_L_plot = {};
M_plot = {};
Re_plot1 = {};
Re_plot2 = {};

k = 1;

for x = C_D_iter
    
C_D = x*ones(1,500)';

    for gamma_0 = gamma_iter

        C_L = zeros(1,500)';
        time = linspace(0, 500, 500)';
        plotting = false;

        disp('Solving Diff Eq'); 
        [M1, V, Z, time, rho, P, T, A, Re] = convergence_numsoln(C_D, C_L, time, gamma_0);

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
        Re = Re(1:clip_index);

        C_D_plot{k} = C_D;
        C_L_plot{k} = C_L;
        M_plot{k} = M1;
        Re_plot{k} = Re;
        disp(k)
        k = k+1;
    end
k = k + 1;
end
disp('done')
%% Plotting
%{
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
%}
figure(1)
set(gcf,'color','w');
for i = 1:8
    plot(M_plot{i}, Re_plot{i});
    hold on; 
end
legend('\gamma = -1.4 C_D = 0.8', '\gamma = -10 C_D = 0.8', '\gamma = -40 C_D = 0.8', '\gamma = -60 C_D = 0.8', '\gamma = -1.4 C_D = 1.4', '\gamma = -10 C_D = 1.4', '\gamma = -40 C_D = 1.4', '\gamma = -60 C_D = 1.4');
ylabel('Reynolds Number');
xlabel('Mach Number');
title('Reynolds Number vs Mach Number');


