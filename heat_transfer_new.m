n = 100; %discretization
R_N = 0.272; %m

thickness = [0.001, 0.001, 0.001, 0.001];
total_thickness = sum(thickness); %m
alphas = [alpha_FW12, alpha_FW12, alpha_FW12, alpha_FW12]; % m^2/s
lambdas = [lambda_FW12, lambda_FW12, lambda_FW12, lambda_FW12]; %W/m/K

eps = 0.9; 
sigma = 5.6695e-8;

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

traj = table2array(readtable('OptimizedPathDatabase.xlsx'));
time = traj(3:end,1);
V = traj(3:end,2);
alpha_r = traj(3:end,3);
q = traj(3:end,4);
X = traj(3:end,5);
Z = traj(3:end,6);
theta_r = traj(3:end,7);
gamma_r = traj(3:end,8);

gamma_atmos = 1.4; 
R_atmos = 287; %J/kg/K


[Z_total, Z_L, Z_U, T, P, rho, c, g_atmos, mu, nu, k, n_eh, n_sum]  = atmo(120, 0.01, 1);

% flip trajectory
Z_total = flip(Z_total)*1000; %convert to m
P = flip(P);
rho = flip(rho);
T = flip(T); 

T_traj = interp1(Z_total, T, Z);
P_traj = interp1(Z_total, P, Z);
rho_traj = interp1(Z_total, rho, Z);

r = 0.71; 
cp = 1.4;

T_aw = T_traj + r*V.^2/(2*cp);

%% Setup integration parameters 

dx = total_thickness/n;
dt = 0.9*0.5*dx^2./max(alphas); 

total_time = time(length(time));
integration_len = round(total_time/dt);
time_integ = (0:integration_len-1)*dt;

T_aw_integ = interp1(time, T_aw, time_integ); 
V_integ = interp1(time, V, time_integ); 
rho_integ = interp1(time, rho_traj, time_integ); 

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
    b_next = 1/2*alphas(i+1)*dt/(dx^2);
    A_total(indices(i),indices(i)+1) = -b;
    A_total(indices(i)+1,indices(i)) = -b;
    A_total(indices(i), indices(i)) = 1 + b + b_next; 
end
%%
n = length(A_total);
T = ones(n,1)*225; %K initial atmospheric temperature at 120 km across traj
T_front = [];
T_back = [];
T_13 = []; 
T_23 = []; 

qw_integ = []; 

if at_stagnation
            N = 0.5; 
            M = 3; 
            qw_integ(1) = 1.83e-4 * R_N^(-1/2) * (1 - T(1)/T_aw_integ(1)) * rho_integ(1)^N * V_integ(1)^M;
else 
            curvilinear_abscissa = 0.56; 
            N = 0.5; 
            M = 3.2;
            qw_integ(idx+1) = 4.03e-5 * cos(pi/4)^(0.5) * sin(pi/4) * curvilinear_abscissa^(-1/2) * (1 - T(1)/T_aw_integ(1)) * rho_integ(1)^N * V_integ(1)^M; 
end

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
    
    for j = 1:length(sizes)-1
        b  = 1/2*alphas(j)*dt/(dx^2);
        b_next = 1/2*alphas(j+1)*dt/(dx^2);
        q(indices(j)) = b*T(i-1)+(1-b-b_next)*T(i)+b_next*T(i+1);
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
            M = 3.2;
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

