function [time, V, gamma_r, Z, T_traj, rho_traj, P_traj, Mach, C_p0, T_w_TM_stag, T_aw, q_w_stag_TM, C_p, T_w_TM, T_aw_TM, q_w_TM] = traj_params_heat(sigma_inp, epsilon)
%Read in variables
traj = table2array(readtable('OptimizedPathDatabase.xlsx'));
time = traj(3:end,1);
V = traj(3:end,2);
alpha_r = traj(3:end,3);
q = traj(3:end,4);
X = traj(3:end,5);
Z = traj(3:end,6);
theta_r = traj(3:end,7);
gamma_r = traj(3:end,8);

% Calculate Mach, Rho, P_inf along traj
% First need T from atmo function
gamma_atmos = 1.4;
R_atmos = 287; %J/kg/K

% sigma = 5.6695e-8;
% eps = 0.9;

sigma = sigma_inp; 
eps = epsilon; 

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



% Find Cp_0 along the Flight Path
gamma = gamma_atmos;
M1 = Mach;
% nose pressure, contingent on shock theory and M1 dependance
P2_P1 = 1 + 2*gamma/(gamma+1)*(M1.^2-1);
Pstag_P2 = ((gamma+1)^2*M1.^2./(4*gamma*M1.^2-2*(gamma-1))).^3;
C_p0 = 2./(gamma*M1.^2).*(P2_P1 .* Pstag_P2-1);

% Calculate wall temperature using TM at stagnation point (re-radiation = convec + rad)
R_N = 0.272; %m
r = 0.71^0.5; %laminar flow
T_aw = T_traj + r*V.^2./(2*C_p0);

A_q = 1.83e-4 * R_N^(-1/2) * rho_traj.^(1/2) .* V.^(3); 
B_q = A_q ./ T_aw; 


T_w_TM_stag = []; 

for i = 1:length(T_aw)
    func_haha = [1 0 0 B_q(i)/(eps*sigma) -1*A_q(i)/(eps*sigma)]; 
    rn = roots(func_haha); 
    temp = rn(rn == real(rn)); 
    T_w_TM_stag(i) = temp(temp > 0);
end



% Tauber Menees Stagnation flux
N = 0.5; 
M = 3; 

C = []; 
q_w_stag_TM = []; 

for i = 1:length(T_aw)
    C(i) = 1.83e-4 * R_N^(-1/2) * (1 - T_w_TM_stag(i)/T_aw(i)); 
    q_w_stag_TM(i) = C(i) * rho_traj(i)^N * V(i)^M; 
end


% Q_w using Tauber Menees in the x_z plane
TR = stlread('CAD_capsule_3.stl');
pts = incenter(TR);
F = faceNormal(TR);


% vect_indices_1 = find(pts(:,2) > -5);
% vect_indices_2 = find(pts(:,2) < 5); 
% vect_indices_3 = find(pts(:,1) > 0); 
% %vect_indices_4 = find(pts(:,3) > 0); 
% [vect_indices, eh] = intersect(vect_indices_1, vect_indices_2);
% [vect_indices, eh] = intersect(vect_indices, vect_indices_3); 
% %[vect_indices, eh] = intersect(vect_indices, vect_indices_4); 

point_index = 7659; 
curvilinear_abscissa = 0.56; 

C_p = []; 
T_w_TM = []; 
T_aw_TM = []; 
q_w_TM = []; 


for i = 1:length(V)
    V_mag = V(i); 
    alpha = alpha_r(i);

    vx = V_mag * sin(alpha);
    vy = 0;
    vz = V_mag * cos(alpha);

    V_vec = [vx, vy, vz];

    A_q = 1.83e-4 * R_N^(-1/2) * rho_traj(i)^(1/2) * V_mag^(3);

    N = F(point_index, :); 
    sin_theta = sum((V_vec.*N)/V_mag);

    C_p(i) = C_p0(i) * sin_theta^2; 
    T_aw_TM(i) = T_traj(i) + r * V_mag^2/(2*C_p(i)); 

    B_q = A_q/T_aw_TM(i); 

    func_haha = [1 0 0 B_q/(eps*sigma) -1*A_q/(eps*sigma)]; 
    rn = roots(func_haha); 
    temp = rn(rn == real(rn)); 
    T_w_TM(i) = temp(temp > 0);



    C = 4.03e-5 * cos(pi/4)^(0.5) * sin(pi/4) * curvilinear_abscissa^(-1/2) * (1 - T_w_TM(i)/T_aw_TM(i)); 
    N = 0.5; 
    M = 3.2; 
    q_w_TM(i) = C * rho_traj(i)^N * V_mag^M; 
    

end

end

