%qw = q_w_stag_TM(1:12000); 
n = 100; %discretization

R_N = 0.272; %m

thickness = [0.001, 0.001, 0.001, 0.001, 0.001, 0.0005, 0.0005];
total_thickness = sum(thickness); %m
alphas = [alpha_FW12, alpha_FW12, alpha_FW12, alpha_FW12, alpha_FW12, alpha_FW12, alpha_FW12]; % m^2/s
lambdas = [lambda_FW12, lambda_FW12, lambda_FW12, lambda_FW12, lambda_FW12, lambda_FW12, lambda_FW12]; %W/m/K

eps = 0.9; 
sigma = 5.6695e-8;
[time, V, gamma_r, Z, T_traj, rho_traj, P_traj, Mach, C_p0, T_w_TM_stag, T_aw_TM_stag, q_w_stag_TM, C_p, T_w_TM, T_aw_TM, q_w_TM] = traj_params_heat(sigma, eps);

% set this to q_w_stag_TM if looking at stagnation point
% set this to q_w_TM if looking at point midway on cone
qw = q_w_stag_TM(1:12000); 
T_w = T_w_TM_stag(1:12000); 
T_aw = T_aw_TM_stag(1:12000); 
at_stagnation = true; 

cumthick = 0;
cum_thickness = ones(7,1);
sizes = zeros(7,1);
indices = zeros(7,1);
current_pos = 0;

for i= 1:length(thickness)
    cumthick = cumthick + thickness(i);
    cum_thickness(i) = cumthick;
    current_pos = thickness(i)/total_thickness*n;
    sizes(i) = round(current_pos);
    indices(i) = sum(sizes);
end

%% Setup integration parameters 

dx = total_thickness/n;
dt = 0.9*0.5*dx^2./max(alphas); 

total_time = time(12000);
integration_len = round(total_time/dt);
time_integ = (0:integration_len-1)*dt;
qw_integ = interp1(time(1:12000), qw, time_integ); 
T_aw_integ = interp1(time(1:12000), T_aw(1:12000), time_integ); 
V_integ = interp1(time(1:12000), V(1:12000), time_integ); 
rho_integ = interp1(time(1:12000), rho_traj(1:12000), time_integ); 
T_w_integ = interp1(time(1:12000), T_w, time_integ); 

%% Setup matrixes

A_total = [];
for i = 1:length(sizes)
    b = ones(sizes(i),1)*1/2*alphas(i)*dt/(dx^2);
    A = spdiags([-b 1+2*b -b],-1:1, sizes(i),sizes(i));
    if i == 1
        A(1,1) = -1;
        A(1,2) = 1;
    elseif i == length(sizes)
        A(sizes(i),sizes(i)-1) = -1;
        A(sizes(i),sizes(i)) = 1;
    end
    A_total = blkdiag(A_total, full(A));
end

for i= 1:length(sizes)-1
    b  = 1/2*alphas(i)*dt/(dx^2);
    A_total(indices(i),indices(i)+1) = -b;
    A_total(indices(i)+1,indices(i)) = -b;
end
%%
n = length(A_total);
T = ones(n,1)*225; %K initial atmospheric temperature at 120 km across traj
T_front = [];
T_back = [];
T_13 = []; 
T_23 = []; 

for idx=1:integration_len
    % set q matrix each time with updated stangation flux, qw(idx)
    q = [];
    for i = 1:n
        if i == 1
            q(i) = -dx/lambdas(1)*(qw_integ(idx)-sigma*eps*T(1)^4);
        elseif i == n
            q(i) = -dx/lambdas(7)*sigma*eps*T(n)^4;
        else
            if i <= indices(1)
                b_mat = 0.5*alphas(1)*dt/(dx^2);
            elseif (indices(1) < i) && (i <= indices(2))
                b_mat = 0.5*alphas(2)*dt/(dx^2);
            elseif (indices(2) < i) && (i <= indices(3))
                b_mat = 0.5*alphas(3)*dt/(dx^2);
            elseif (indices(3) < i) && (i <= indices(4))
                b_mat = 0.5*alphas(4)*dt/(dx^2);
            elseif (indices(4) < i) && (i <= indices(5))
                b_mat = 0.5*alphas(5)*dt/(dx^2);
            elseif (indices(5) < i) && (i <= indices(6))
                b_mat = 0.5*alphas(6)*dt/(dx^2);
            else
                b_mat = 0.5*alphas(7)*dt/(dx^2);
            end
            q(i) = b_mat*T(i-1)+(1-2*b_mat)*T(i)+b_mat*T(i+1);
        end
    end
    
    T = A_total\transpose(q);
    T_front(idx) = T(1);
    T_13(idx) = T(round(n/3)); 
    T_23(idx) = T(round(2*n/3)); 
    T_back(idx) = T(n);

    if at_stagnation
        if idx ~= integration_len
            disp('here')
            N = 0.5; 
            M = 3; 
            qw_integ(idx + 1) = 1.83e-4 * R_N^(-1/2) * (1 - T_front(idx)/T_aw_integ(idx)) * rho_integ(idx)^N * V_integ(idx)^M;
        end
    else 
        if idx~= integration_len
            curvilinear_abscissa = 0.56; 
            N = 0.5; 
            M = 3.2
            qw_integ(idx+1) = 4.03e-5 * cos(pi/4)^(0.5) * sin(pi/4) * curvilinear_abscissa^(-1/2) * (1 - T_front(idx)/T_aw_integ(idx)) * rho_integ(idx)^N * V_integ(idx)^M; 
        end
    end

end

%% Plotting

baseline = ones(1,length(T_front))*(70+273.15); %back wall
figure(1)
plot(T_front)
hold on
plot(T_back)
hold on
plot(T_13)
hold on
plot(T_23)
hold on
plot(baseline, '--')
legend('front', '6 mm','2 mm','4 mm','70 C')
set(gcf,'color','w');

%%
figure(2)
plot(q_w_TM(1:12000))
hold on
plot(q_w_stag_TM(1:12000))
legend('MidPoint', 'Cone')
set(gcf,'color','w');


