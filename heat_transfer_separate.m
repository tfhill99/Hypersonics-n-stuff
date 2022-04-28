n = 50; %discretization per layer
R_N = 0.272; %m

thickness = [0.002, 0.002, 0.02, 0.02]; %m
total_thickness = sum(thickness); %m
alphas = [alpha_Rescor311, alpha_Rescor310M, alpha_Rescor310M, alpha_Rescor310M]; % m^2/s
lambdas = [lambda_Rescor311, lambda_Rescor310M, lambda_Rescor310M, lambda_Rescor310M]; %W/m/K

eps = 0.9; 
sigma = 5.6695e-8;

cumthick = 0;
cum_thickness = zeros(4,1);
sizes = zeros(4,1);
indices = zeros(4,1);
current_pos = 0;

for i= 1:length(thickness)
    cumthick = cumthick + thickness(i);
    cum_thickness(i) = cumthick;
    current_pos = thickness(i)/total_thickness*n;
    sizes(i) = round(current_pos);
    indices(i) = sum(sizes);
end

positions = [1/3 2/3]; 
indices = []; 
for i = 1:length(positions)
    curr_thickness = positions(i) * total_thickness; 
    b = cum_thickness(1:(find(cum_thickness > curr_thickness)-1)); 
    excessive_thickness = curr_thickness - cum_thickness(length(b)); 
    indices(i) = length(b) * n + round(n/thickness(length(b) + 1) * excessive_thickness);
end

traj = table2array(readtable('OptimizedPathDatabase.xlsx'));
time = traj(3:end,1);
V = traj(3:end,2);
alpha_r = traj(3:end,3);
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

dx = thickness/n;
dt = 0.9*0.5*(dx.^2)./alphas; 

total_time = time(length(time));
integration_len = round(total_time/min(dt));
time_integ = (0:integration_len-1)*min(dt);

T_aw_integ = interp1(time, T_aw, time_integ); 
V_integ = interp1(time, V, time_integ); 
rho_integ = interp1(time, rho_traj, time_integ); 

%% Setup matrixes
As = {}; 

for j = 1:length(thickness)
        b = ones(n,1)*1/2*alphas(j)*dt(j)/(dx(j)^2);
        A = spdiags([-b 1+2*b -b],-1:1, n, n);
            A(1,1) = -1;
            A(1,2) = 1;
            A(n, n-1) = -1;
            A(n, n) = 1;
        As{j} = full(A);
end
% 
% for i= 1:length(sizes)-1
%     b  = 1/2*alphas(i)*dt/(dx^2);
%     b_next = 1/2*alphas(i+1)*dt/(dx^2);
%     A_total(indices(i),indices(i)+1) = -b;
%     A_total(indices(i)+1,indices(i)) = -b;
%     A_total(indices(i), indices(i)) = 1 + b + b_next; 
% end
%%
%n = length(A_total);
T = ones(4*n,1)*225; %K initial atmospheric temperature at 120 km across traj
T_front = [];
T_back = [];
T_blayer1 = [];
T_blayer2 = [];
T_blayer3 = [];
at_stagnation = true;

qw_integ = []; 

if at_stagnation
    N = 0.5; 
    M = 3; 
    qw_integ(1) = 1.83e-4 * R_N^(-1/2) * (1 - T(1)/T_aw_integ(1)) * rho_integ(1)^N * V_integ(1)^M;
else 
    curvilinear_abscissa = 0.56; 
    N = 0.5; 
    M = 3.2;
    qw_integ(1) = 4.03e-5 * cos(pi/4)^(0.5) * sin(pi/4) * curvilinear_abscissa^(-1/2) * (1 - T(1)/T_aw_integ(1)) * rho_integ(1)^N * V_integ(1)^M; 
end

for idx=1:integration_len
    % set q matrix each time with updated stangation flux, qw(idx)
    q = [];
    for i = 1:n
        if i == 1
            q(i) = -dx(1)/lambdas(1)*(qw_integ(idx)-sigma*eps*T(1)^4);
        elseif i == n
            q(i) = -dx(1)/lambdas(1)*sigma*eps*T(n)^4;
        else
            b_mat = 0.5*alphas(1)*dt(1)/(dx(1)^2);
            q(i) = b_mat*T(i-1)+(1-2*b_mat)*T(i)+b_mat*T(i+1);
        end
    end

    b = 0.5 * alphas(1) * dt(1)/dx(1)^2; 
    b = 0.5 * alphas(2) * dt(2)/dx(2)^2; 
    q(n) = b*T(n-1) + (1-b-b_next)*T(n) + b_next*T(n+1);

    for i = n+1:2*n
        if i == n+1
            q(i) = -dx(2)/lambdas(2)*(q(n)-sigma*eps*T(n+1)^4);
        elseif i == 2*n
            q(i) = -dx(2)/lambdas(2)*sigma*eps*T(2*n)^4;
        else
            b_mat = 0.5*alphas(2)*dt(2)/(dx(2)^2);
            q(i) = b_mat*T(i-1)+(1-2*b_mat)*T(i)+b_mat*T(i+1);
        end
    end

    b = 0.5 * alphas(2) * dt(2)/dx(2)^2; 
    b = 0.5 * alphas(3) * dt(3)/dx(3)^2; 
    q(2*n) = b*T(2*n-1) + (1-b-b_next)*T(2*n) + b_next*T(2*n+1);

    for i = 2*n+1:3*n
        if i == 2*n + 1
            q(i) = -dx(3)/lambdas(3)*(q(2*n)-sigma*eps*T(2*n+1)^4);
        elseif i == 3*n
            q(i) = -dx(3)/lambdas(3)*sigma*eps*T(3*n)^4;
        else
            b_mat = 0.5*alphas(3)*dt(3)/(dx(3)^2);
            q(i) = b_mat*T(i-1)+(1-2*b_mat)*T(i)+b_mat*T(i+1);
        end
    end

    b = 0.5 * alphas(3) * dt(3)/dx(3)^2; 
    b = 0.5 * alphas(4) * dt(4)/dx(4)^2; 
    q(3*n) = b*T(3*n-1) + (1-b-b_next)*T(3*n) + b_next*T(3*n+1);

    for i = 3*n+1:4*n
        if i == 3*n+1
            q(i) = -dx(4)/lambdas(4)*(q(3*n)-sigma*eps*T(3*n+1)^4);
        elseif i == 4*n
            q(i) = -dx(4)/lambdas(4)*sigma*eps*T(4*n)^4;
        else
            b_mat = 0.5*alphas(4)*dt(4)/(dx(4)^2);
            q(i) = b_mat*T(i-1)+(1-2*b_mat)*T(i)+b_mat*T(i+1);
        end
    end    
    

    T1 = As{1}\transpose(q(1:n));
    T2 = As{2}\transpose(q((n+1):2*n));
    T3 = As{3}\transpose(q((2*n+1):3*n));
    T4 = As{4}\transpose(q((3*n+1):4*n));
    T = [T1; T2; T3; T4]; 

    T_front(idx) = T(1);
    T_blayer1(idx) = T(indices(1)); 
    T_blayer2(idx) = T(indices(2)); 
    %T_blayer3(idx) = T(indices(3)); 
    T_back(idx) = T(4*n);

    if at_stagnation
        if idx ~= integration_len
            N = 0.5; 
            M = 3; 
            qw_integ(idx + 1) = 1.83e-4 * R_N^(-1/2) * (1 - T_front(idx)/T_aw_integ(idx)) * rho_integ(idx)^N * V_integ(idx)^M;
        end
    else 
        if idx ~= integration_len
            curvilinear_abscissa = 0.56; 
            N = 0.5; 
            M = 3.2;
            qw_integ(idx + 1) = 4.03e-5 * cos(pi/4)^(0.5) * sin(pi/4) * curvilinear_abscissa^(-1/2) * (1 - T_front(idx)/T_aw_integ(idx)) * rho_integ(idx)^N * V_integ(idx)^M; 
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
plot(T_blayer1)
hold on
plot(T_blayer2)
hold on
%plot(T_blayer3)
%hold on
plot(baseline, '--')
legend('front', 'back', 'layer 1 back','layer 2 back', '70 C')
set(gcf,'color','w');

%%
figure(2)
plot(time_integ, qw_integ);

