%qw = q_w_stag_TM(1:12000); 
n = 100; %discretization
%total_time = time(12000); %steps of solving

thickness = [0.001, 0.001, 0.001, 0.001, 0.001, 0.0005, 0.0005];
total_thickness = sum(thickness); %m
alphas = [alpha_FW12, alpha_Rescor310M, alpha_Rescor311, alpha_nextel, alpha_sigratherm, alpha_pyrogel, alpha_pyrogel]; % m^2/s
lambdas = [lambda_FW12, lambda_Rescor310M, lambda_Rescor311, lambda_nextel, lambda_sigratherm, lambda_pyrogel, lambda_pyrogel]; %W/m/K

eps = 0.9; 
sigma = 5.6695e-8;
[time, V, gamma_r, Z, T_traj, rho_traj, P_traj, Mach, C_p0, T_w_TM_stag, q_w_stag_TM, C_p, T_w_TM, q_w_TM] = traj_params_heat(sigma, eps);

qw = q_w_stag_TM; 

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

% sigma = 5.6695e-8; 
% eps = 0.9; 

dx = total_thickness/n;
dts = 0.7*0.5*dx^2./alphas; %5e-6

%% Setup matrixes

A_total = [];
for i = 1:length(sizes)
    b = ones(sizes(i),1)*1/2*alphas(i)*dts(i)/(dx^2);
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
    b  = 1/2*alphas(i)*dts(i)/(dx^2);
    A_total(indices(i),indices(i)+1) = -b;
    A_total(indices(i)+1,indices(i)) = -b;
end
%%

T = ones(n,1)*225; %K initial atmospheric temperature at 120 km across traj
T_front = [];
T_back = [];
T_13 = []; 
T_23 = []; 
n = length(A_total);

for idx=1:12000
    % set q matrix each time with updated stangation flux, qw(idx)
    q = [];
    for i = 1:n
        if i == 1
            q(i) = -dx/lambdas(1)*(qw(idx)-sigma*eps*T(1)^4);
        elseif i == n
            q(i) = -dx/lambdas(7)*sigma*eps*T(n)^4;
        else
            if i <= indices(1)
                b_mat = 0.5*alphas(1)*dts(1)/(dx^2);
            elseif (indices(1) < i) && (i <= indices(2))
                b_mat = 0.5*alphas(2)*dts(2)/(dx^2);
            elseif (indices(2) < i) && (i <= indices(3))
                b_mat = 0.5*alphas(3)*dts(3)/(dx^2);
            elseif (indices(3) < i) && (i <= indices(4))
                b_mat = 0.5*alphas(4)*dts(4)/(dx^2);
            elseif (indices(4) < i) && (i <= indices(5))
                b_mat = 0.5*alphas(5)*dts(5)/(dx^2);
            elseif (indices(5) < i) && (i <= indices(6))
                b_mat = 0.5*alphas(6)*dts(6)/(dx^2);
            else
                b_mat = 0.5*alphas(7)*dts(7)/(dx^2);
            end
            disp(i)
            q(i) = b_mat*T(i-1)+(1-2*b_mat)*T(i)+b_mat*T(i+1);
        end
    end
    
    T = A_total\transpose(q);
    T_front(idx) = T(1);
    T_13(idx) = T(round(n/3)); 
    T_23(idx) = T(round(2*n/3)); 
    T_back(idx) = T(n);
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
legend('front', 'back','1/3','2/3','limit')
set(gcf,'color','w');


