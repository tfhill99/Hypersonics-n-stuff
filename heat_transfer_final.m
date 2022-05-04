function[temperatures] = heat_transfer_new(thickness)

n = 300; %discretization  % 100-400 for rigid, 300 for flexible
R_N = 0.272; %m

%% Rigid TPS Materials
[alpha_FW12, alpha_Rescor310M, ~, alpha_Intek1120, ~, ~, ~, lambda_FW12, lambda_Rescor310M, ~, lambda_Intek1120, ~, ~, ~] = TPS_materials();
%% Flexible TPS Materials
[~, ~, ~ , ~, alpha_nextel, alpha_pyrogel, alpha_sigratherm, ~, ~, ~, ~, lambda_nextel, lambda_pyrogel, lambda_sigratherm] = TPS_materials();

%%
%thickness = [0.003, 0.003, 0.04, 0.002]; %thickness in m
total_thickness = sum(thickness); %m
alphas = [alpha_FW12, alpha_Rescor310M, alpha_Intek1120, alpha_FW12]; %m^2/s RIGID
lambdas = [lambda_FW12, lambda_Rescor310M, lambda_Intek1120, lambda_FW12];%W/m/K RIGID
%alphas = [alpha_nextel, alpha_pyrogel, alpha_sigratherm, alpha_nextel]; % FLEXIBLE
%lambdas = [lambda_nextel, lambda_pyrogel, lambda_sigratherm, lambda_nextel]; % FLEXIBLE

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
cumthick %display new cumthick each iteration

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

[Z_total, Z_L, Z_U, Temp, P, rho, c, g_atmos, mu, nu, k, n_eh, n_sum]  = atmo(120, 0.01, 1);

% flip trajectory
Z_total = flip(Z_total)*1000; %convert to m
rho = flip(rho);
Temp = flip(Temp); 

T_traj = interp1(Z_total, Temp, Z);
rho_traj = interp1(Z_total, rho, Z);

r = 0.71; 
cp = 1.4;

T_aw = T_traj + r*V.^2/(2*cp);

%% Setup integration parameters 

dx = total_thickness/n;
dt = 0.999*0.5*dx^2/max(alphas); 

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
    %b_avg = 0.5*(b+b_next);
    %A_total(indices(i)-1,indices(i)) = -b_avg;
    %A_total(indices(i),indices(i)-1) = -b_avg;
    A_total(indices(i),indices(i)+1) = -b_next;
    A_total(indices(i)+1,indices(i)) = -b_next;
    A_total(indices(i), indices(i)) = 1 + b + b_next; 
end

%%
n = length(A_total);
%T = ones(n,1); %K initial atmospheric temperature at 120 km across traj
%T = T_traj;
T = [n, 1];
T(1) = T_traj(1);
for i = 2:n-1
    T(i) = T_traj(fix(length(T_traj)/n)*i);
end
T(n) = T_traj(n);

T_front = [];
T_back = [];
T_blayer1 = [];
T_blayer2 = [];
T_blayer3 = [];
T_thickness = {};

at_stagnation = true; % toggle for rigid=true, flexible=false

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
            q(i) = -dx/lambdas(1)*(qw_integ(idx)-sigma*eps*T(1)^4);
        elseif i == n
            q(i) = -dx/lambdas(4)*sigma*eps*T(n)^4;
        else
            if i <= indices(1)
                b_mat = 0.5*alphas(1)*dt/(dx^2);
            elseif (indices(1) < i) && (i <= indices(2))
                b_mat = 0.5*alphas(2)*dt/(dx^2);
            elseif (indices(2) < i) && (i <= indices(3))
                b_mat = 0.5*alphas(3)*dt/(dx^2);
            else
                b_mat = 0.5*alphas(4)*dt/(dx^2);
            end
            q(i) = b_mat*T(i-1)+(1-2*b_mat)*T(i)+b_mat*T(i+1);
        end
    end
    
    for j = 1:length(sizes)-1
        b  = 1/2*alphas(j)*dt/(dx^2);
        b_next = 1/2*alphas(j+1)*dt/(dx^2);
        q(indices(j)) = b*T(indices(j)-1)+(1-b-b_next)*T(indices(j))+b_next*T(indices(j)+1);
    end

    T = A_total\transpose(q);
    if mod(idx,100) == 0
        T_thickness{idx/100} = T;
    end
    T_front(idx) = T(1);;
    T_blayer1(idx) = T(indices(1)); 
    T_blayer2(idx) = T(round((indices(2)-indices(1))/2) + indices(1)); 
    T_blayer3(idx) = T(round((indices(3)-indices(2))/2) + indices(2)); 
    T_back(idx) = T(n);

    if at_stagnation
        if idx ~= integration_len
            N = 0.5; 
            M = 3; 
            qw_integ(idx + 1) = 1.83e-4 * R_N^(-1/2) * (1 - (T_front(idx))/T_aw_integ(idx)) * rho_integ(idx)^N * V_integ(idx)^M;
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

temperatures = [max(T_blayer1); max(T_blayer2); max(T_blayer3); max(T_back)];

%% Plotting along time

baseline = ones(1,length(T_front))*(70+273.15); %back wall
opT = ones(1,length(T_front))*(2000); %front wall
figure(1)
plot(time_integ,T_front)
hold on
plot(time_integ,T_blayer1)
hold on
plot(time_integ,T_blayer2)
hold on
plot(time_integ,T_blayer3)
hold on
plot(time_integ,T_back)
hold on
plot(opT, '--','LineWidth',1)
hold on
plot(baseline, '--','LineWidth',1)
title('TPS Temperature over Flight Path')
ylim([200 2100])
xlim([0 time_integ(end)])
ylabel('Temperature (K)')
xlabel('Time (s)')
legend('T_{front wall}', 'T_{layer 1} (FW12)','T_{layer 2} (FW12 + Rescor310M)','T_{layer 3} (FW12 + Rescor310M + Intek1120)', 'T_{back wall} (FW12 + Rescor310M + Intek1120 + FW12)','T_{operational}', 'T_{payload}', 'FontSize', 7)
set(gcf,'color','w');

%% Plotting along thickness
x = linspace(0, cumthick, n)*1000;
baseline_t = ones(1,length(x))*(70+273.15); %back wall
figure(2)
plot(x, T_thickness{1})
hold on
plot(x, T_thickness{10})
hold on
plot(x, T_thickness{30})
hold on
plot(x, T_thickness{50})
hold on
plot(x, T_thickness{100})
hold on
plot(x, T_thickness{length(T_thickness)})
hold on
plot(x, baseline_t, '--', 'LineWidth',1)
title('TPS Temperature across Thickness')
ylabel('Temperature (K)')
xlabel('Distance (mm)')
xlim([0 cumthick*1000])
legend('t_{0 s}', 't_{41.1 s}','t_{123 s}','t_{205.5 s}', 't_{411.5 s}', 't_{1012 s}','T_{payload}')
set(gcf,'color','w');

%% plotting materials vs thickness
x = linspace(0, cumthick, n)*1000;
baseline_t = ones(1,length(x))*(70+273.15); %back wall
figure(3)
plot(0, max(T_front),'*')
hold on
plot(thickness(1)*1000, max(T_blayer1), '*')
plot((thickness(1) + thickness(2))*1000, max(T_blayer2), '*')
plot((thickness(1) + thickness(2) + thickness(3))*1000, max(T_blayer3), '*')
plot((thickness(1) + thickness(2) + thickness(3) + thickness(4))*1000, max(T_back), '*')
title('TPS Temperature across Thickness')
ylabel('Temperature (K)')
xlabel('Distance (mm)')
xlim([0 cumthick*1000])
legend('t_{front}', 't_{1}','t_{2}','t_{3}', 't_{back}')
set(gcf,'color','w');


end
