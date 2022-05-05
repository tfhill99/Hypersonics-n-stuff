clear; 
clc; 
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
T = flip(T); 

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
plot(Mach, Z/1000)
hold on
plot(Load_Factor, Z_acc/1000)
title('Trajectory Properties v/s Altitude (Optimized Flight Path')
legend('Mach', 'Load Factor')
xlabel('Property')
ylabel('Altitude (km)')

%% Find Cp_0 along the Flight Path
gamma = gamma_atmos;
M1 = Mach;
% nose pressure, contingent on shock theory and M1 dependance
P2_P1 = 1 + 2*gamma/(gamma+1)*(M1.^2-1);
Pstag_P2 = ((gamma+1)^2*M1.^2./(4*gamma*M1.^2-2*(gamma-1))).^3;
C_p0 = 2./(gamma*M1.^2).*(P2_P1 .* Pstag_P2-1);

%% Calculate wall temperature using TM at stagnation point (re-radiation = convec + rad)
R_N = 0.272; %m
r = 0.71^0.5; %laminar flow
R_H_NS_factor = ((gamma_atmos + 1) * Mach.^2)./(2 + (gamma_atmos - 1)*Mach.^2); 
V_NS = V./R_H_NS_factor; 

const_1 = (2 * gamma_atmos)/(gamma_atmos + 1); 
T_RH = T_traj .* (1 + const_1 * (Mach.^2 - 1)) ./ R_H_NS_factor; 

T_aw = T_RH + r*V_NS.^2/(2*1.4);
T_aw_2 = T_traj + r*V.^2/(2*1.4); 

A_q = 1.83e-4 * R_N^(-1/2) * rho_traj.^(1/2) .* V.^(3); 
B_q = A_q ./ T_aw; 
B_q_2 = A_q ./ T_aw_2; 

sigma = 5.6695e-8; 
eps = 0.9; 

T_w_TM = []; 
T_w_TM_2 = []; 

for i = 1:length(T_aw)
    func_haha = [1 0 0 B_q(i)/(eps*sigma) -1*A_q(i)/(eps*sigma)]; 
    rn = roots(func_haha); 
    temp = rn(rn == real(rn)); 
    T_w_TM(i) = temp(temp > 0);

    func_haha = [1 0 0 B_q_2(i)/(eps*sigma) -1*A_q(i)/(eps*sigma)]; 
    rn = roots(func_haha); 
    temp = rn(rn == real(rn)); 
    T_w_TM_2(i) = temp(temp > 0);
end


figure(1)
plot(T_w_TM, Z/1000);
hold on
plot(T_w_TM_2, Z/1000); 
xlabel('T_w (K)'); 
ylabel('Altitude(km)');
title('Tauber Menees Wall Temperature at Stagnation Point');

%% Tauber Menees Stagnation
N = 0.5; 
M = 3; 

C = []; 
q_w_stag_TM = []; 
q_w_stag_TM_2 = []; 

for i = 1:length(T_aw)
    C(i) = 1.83e-4 * R_N^(-1/2) * (1 - T_w_TM(i)/T_aw(i)); 
    q_w_stag_TM(i) = C(i) * rho_traj(i)^N * V(i)^M; 

    C2 = 1.83e-4 * R_N^(-1/2) * (1 - T_w_TM_2(i)/T_aw_2(i)); 
    q_w_stag_TM_2(i) = C2 * rho_traj(i)^N * V(i)^M; 
end

figure(2)
plot(q_w_stag_TM(1:12000)/1000, Z(1:12000)); 
xlabel('q_w (kW/m^2)'); 
ylabel('Altitude(m)');
title('Tauber Menees q_w at Stagnation');

figure(3)
plot(time(1:13000), q_w_stag_TM(1:13000)/1000); 
%hold on
%plot(time(1:12000), q_w_stag_TM_2(1:12000)/1000);
ylabel('q_w (kW/m^2)'); 
xlabel('Time (s)');
%title('Time History of Stagnation q_w Using Tauber Menees');

%% discretization
z_needed = [120 110 100 90 82 78 70 65 62 58 56 52 48 42 35 30] * 1000; 
V = interp1(Z, V, z_needed); 
rho_traj = interp1(Z, rho_traj, z_needed); 
alpha_r = interp1(Z, alpha_r, z_needed); 
q_w_stag_TM_new = interp1(Z, q_w_stag_TM, z_needed); 
time_new = interp1(Z, time, z_needed); 
T_traj = interp1(Z, T_traj, z_needed); 
Mach_new = interp1(Z, Mach, z_needed); 
Z_new = interp1(Z, Z, z_needed); 


%% checking discretization
figure(2)
plot(q_w_stag_TM(1:12000)/1000, Z(1:12000)); 
hold on; 
scatter(q_w_stag_TM_new/1000, Z_new)
xlabel('q_w (kW/m^2)'); 
ylabel('Altitude(m)');
title('Tauber Menees q_w at Stagnation');


figure(3)
plot(time(1:12000), q_w_stag_TM(1:12000)/1000); 
hold on; 
scatter(time_new, q_w_stag_TM_new/1000); 
ylabel('q_w (kW/m^2)'); 
xlabel('Time (s)');
title('Time History of Stagnation q_w Using Tauber Menees');
%% Q_w using Tauber Menees in the x_z plane
TR = stlread('CAD_capsule_3.stl');
pts = incenter(TR);
F = faceNormal(TR);

broad_local_alpha = {}; 
broad_local_C_p = {}; 
broad_local_T_W_TM = {}; 
broad_local_T_aw_TM = {}; 
broad_local_qw = {}; 
broad_local_x = {}; 
broad_local_x_dist = {}; 
broad_local_x_dist_plot = {}; 

vect_indices_1 = find(pts(:,2) > -5);
vect_indices_2 = find(pts(:,2) < 5); 
vect_indices_3 = find(pts(:,1) > 0); 
%vect_indices_4 = find(pts(:,3) > 0); 
[vect_indices, eh] = intersect(vect_indices_1, vect_indices_2);
[vect_indices, eh] = intersect(vect_indices, vect_indices_3); 
%[vect_indices, eh] = intersect(vect_indices, vect_indices_4); 


for i = 1:length(V)
    V_mag = V(i); 
    alpha = alpha_r(i);

    vx = V_mag * sin(alpha);
    vy = 0;
    vz = V_mag * cos(alpha);

    V_vec = [vx, vy, vz];

    A_q = 1.83e-4 * R_N^(-1/2) * rho_traj(i)^(1/2) * V_mag^(3);

    sigma = 5.6695e-8;
    eps = 0.9;

    local_alphas = []; 
    C_p = [];
    T_w = []; 
    T_aw = []; 

    disp(i)

    T_aw_test = []; 
    T_w_test = []; 

    F_val = []; 

    F_test =[]; 

    for j = 1:length(F)
        N = F(j, :); 
        sin_theta = sum((V_vec.*N)/V_mag); 
        local_alphas(j) = asin(sin_theta); 
        if local_alphas(j) <= 0
            C_p(j) = -2/(1.4 * Mach_new(i)^2); 
        else
            C_p(j) = C_p0(i) * sin_theta^2; 
        end

        RH_factor = ((gamma_atmos + 1) * Mach_new(i)^2)/(2 + (gamma_atmos - 1)*Mach(i)^2);

%         beta = obliquerelations('mach', Mach_new(i), 'theta', pi/4, 1.4); 
%         disp(rad2deg(beta))

        beta = deg2rad(64.2); 

        v_n = V_mag * sin(beta); 
        v_t = V_mag * cos(beta); 

        v_n2 = v_n/RH_factor; 

        V_mag_new = norm([v_t, v_n2]); 

        T_RH = T_traj(i) * (1 + const_1 * (Mach_new(i)^2 - 1)) / RH_factor;


        if i < 9
            T_aw(j) = T_RH + r*V_mag_new^2/(2*1.4);
            T_aw_test(j) = T_traj(i) + r*V_mag^2/(2*1.4);
        else
            T_aw(j) = T_traj(i) + r*V_mag^2/(2*1.4);
        end

        B_q = A_q / T_aw(j); 

        func_haha = [1 0 0 B_q/(eps*sigma) -1*A_q/(eps*sigma)]; 
        rn = roots(func_haha); 
        temp = rn(rn == real(rn)); 
        T_w(j) = temp(temp > 0);

        F_val(j) = 1 - T_w(j)/T_aw(j); 

%         B_q = A_q / T_aw_test(j); 
% 
%         func_haha = [1 0 0 B_q/(eps*sigma) -1*A_q/(eps*sigma)]; 
%         rn = roots(func_haha); 
%         temp = rn(rn == real(rn)); 
%         T_w_test(j) = temp(temp > 0);
% 
%         F_test(j) = 1 - T_w_test(j)/T_aw_test(j); 


    end

    disp(V_mag_new)

    broad_local_alpha{i} = local_alphas; 
    broad_local_C_p{i} = C_p; 
    broad_local_T_W_TM{i} = T_w; 
    broad_local_T_aw_TM{i} = T_aw; 

    cross_sec_cp = []; 
    cross_sec_pts = [];
    cross_sec_Tw = []; 
    cross_sec_Taw = []; 
    cross_sec_alpha = []; 

    enter_index = 1; 

    for lc = 1:length(vect_indices)
        k = vect_indices(lc);
        cross_sec_pts(enter_index, :) = pts(k, :); 
        cross_sec_cp(enter_index) = C_p(k); 
        cross_sec_Tw(enter_index) = T_w(k); 
        cross_sec_Taw(enter_index) = T_aw(k); 
        cross_sec_alpha(enter_index) = local_alphas(k); 
        enter_index = enter_index + 1; 
    end

    %cross_sec_cp = reshape(cross_sec_cp, [length(cross_sec_cp), 1]);
    %cross_sec_pts = reshape(cross_sec_pts, [length(cross_sec_pts), 1]);

    %{
    [cross_sec_pts_x_sorted, I] = sort(cross_sec_pts(:,1));
    cspy = cross_sec_pts(:,2); 
    cspz = cross_sec_pts(:,3); 

    cross_sec_pts = [cross_sec_pts_x_sorted cspy(I) cspz(I)]; 
    cross_sec_cp = cross_sec_cp(I); 
    cross_sec_Tw = cross_sec_Tw(I);
    cross_sec_Taw = cross_sec_Taw(I);
    cross_sec_alpha = cross_sec_alpha(I);
    %}

    [cross_sec_pts_z_sorted, I] = sort(cross_sec_pts(:,3));
    cspy = cross_sec_pts(:,2); 
    cspx = cross_sec_pts(:,1); 

    cross_sec_pts = [flip(cspx(I)) flip(cspy(I)) flip(cross_sec_pts_z_sorted)]; 
    cross_sec_cp = flip(cross_sec_cp(I)); 
    cross_sec_Tw = flip(cross_sec_Tw(I));
    cross_sec_Taw = flip(cross_sec_Taw(I));
    cross_sec_alpha = flip(cross_sec_alpha(I));

    x_dist = [0]; 

    cross_sec_qw = [0]; 

    x_dist_plot = [0]; 

    a = 0.855; 
    b = 0.145; 

    last_index = -1; 
    last_point = [];

    starter = 0; 

    for lc2 = 2:length(cross_sec_cp)
        p2 = cross_sec_pts(lc2, :); 
        p1 = cross_sec_pts(lc2-1, :); 
        % x_dist(lc2) = (x_dist(lc2-1)) + (norm(p2-p1))/1000; 
        % x_dist(lc2) = x_dist(lc2) / 1000;  

        phi = abs(atan(abs(p2(1))/(p2(3)-291)));
        
       
        if phi < pi/4 
            cross_sec_qw(lc2) = q_w_stag_TM_new(i) * (a * cos(phi)^1.5 + b); 
            x_dist_plot(lc2) = phi/(2*pi) * (2 * pi * 0.272); 
            last_index = lc2; 
            last_point = p2; 
            x_dist(lc2) = x_dist_plot(lc2); 
        else
                x_dist(lc2) = x_dist_plot(lc2-1) + (norm(p2 - p1, 2))/1000; 
                x_dist_plot(lc2) = x_dist(lc2); 
            
            %x_dist(lc2) = x_dist(lc2-1) + ((p2(1) - p1(1))*sin(pi/4))/1000; 
            %C = 2.2e-5 * cos(cross_sec_alpha(lc2))^2.08 * sin(cross_sec_alpha(lc2))^1.6 * x_dist(lc2)^(-1/5) * (1 - (1.11 * cross_sec_Tw(lc2))/(cross_sec_Taw(lc2)));
            
            % C = 2.2e-5 * cos(pi/4)^2.08 * sin(pi/4)^1.6 * x_dist(lc2)^(-1/5) * (1 - (1.11 * cross_sec_Tw(lc2))/(cross_sec_Taw(lc2)));

            C = 4.03e-5 * cos(pi/4)^(0.5) * sin(pi/4)^(0.5) * x_dist(lc2)^(-0.5) * (1 - cross_sec_Tw(lc2)/cross_sec_Taw(lc2)); 
            N = 0.5; 
            M = 3.2; 
            cross_sec_qw(lc2) = C * rho_traj(i)^N * V_mag^M; 
        end
       
        %C = 2.2e-5 * cos(cross_sec_alpha(lc2))^2.08 * sin(cross_sec_alpha(lc2))^1.6 * x_dist(lc2)^(-1/5) * (1 - (1.11 * cross_sec_Tw(lc2))/(cross_sec_Taw(lc2)));
        %cross_sec_qw(lc2) = C * rho_traj(i)^N * V_mag^M;

        
    end

    broad_local_qw{i} = cross_sec_qw; 
    broad_local_x_dist{i} = x_dist; 
    broad_local_x{i} = cross_sec_pts(:, 1); 
    broad_local_x_dist_plot{i} = x_dist_plot; 
    
end

%%
qw_curr = broad_local_qw{7}; 
x_curr = broad_local_x{7}; 
x_dist_curr = broad_local_x_dist_plot{7};

figure(10); 
scatter(x_dist_curr, qw_curr/1000);
xlabel('Curvilinear Abscissa (m)'); 
ylabel('q_w (kW/m^2)'); 
title('q_w across cross section at 70 km Altitude');

figure(11)
scatter(broad_local_x_dist_plot{3}, broad_local_qw{3}/1000); 
hold on; 
scatter(broad_local_x_dist_plot{12}, broad_local_qw{12}/1000); 
xlabel('Curvilinear Abscissa (m)'); 
ylabel('q_w (kW/m^2)'); 
title('q_w across cross section at Different Altitudes');
legend('100 km', '52 km'); 
%%
T_w_curr = broad_local_T_W_TM{7}; 

figure(1); 
trisurf(TR, T_w_curr, 'EdgeColor', 'none'); 
title('T_w Distribution at 70 km altitude'); 
hold on;
axis equal; 
colorbar; 
colormapeditor; 


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
