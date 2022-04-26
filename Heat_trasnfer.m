total_thickness = 0.0060; %m
qw = q_w_stag_TM(1:12000); 
n = 200; %discretization
%total_time = time(12000); %steps of solving

thickness = [0.001, 0.001, 0.001, 0.001, 0.001, 0.0005, 0.0005];
alphas = [alpha_FW12, alpha_FW12, alpha_FW12, alpha_FW12, alpha_FW12, alpha_FW12, alpha_FW12]; % m^2/s
lambdas = [lambda_FW12, lambda_FW12, lambda_FW12, lambda_FW12, lambda_FW12, lambda_FW12, lambda_FW12]; %W/m/K

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

sigma = 5.6695e-8; 
eps = 0.9; 

%dxs = thickness/n;
dx = total_thickness/n;
dts = (0.9*0.5)*dx^2./alphas;

%% Setup matrixes

A_total = [];
for i = 1:length(sizes)
    b = ones(sizes(i),1)*1/2*alphas(i)*dts(i)/(dx^2);
    A = spdiags([-b 1+2*b -b],-1:1,sizes(i),sizes(i));
    if i == 1
        A(1,1) = -1;
        A(1,2) = 1;
    elseif i==length(sizes)
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
            q(i) = b_mat*T(i-1)+(1-2*b_mat)*T(i)+b_mat*T(i+1);
        end
    end
    
    T = A_total\transpose(q);
    T_front(idx) = T(1);
    T_13(idx) = T(int16(n/3)); 
    T_23(idx) = T(int16(2*n/3)); 
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
plot(baseline)
legend('front', 'back','1/3','2/3','limit')
set(gcf,'color','w');


