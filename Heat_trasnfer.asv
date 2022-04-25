thickness = 0.060; %m
qw = q_w_stag_TM(1:12000); 
steps = 200; %discretization
total_time = time(12000); %steps of solving

dx = thickness/steps;
dt = (0.9*0.5)*dx^2/alpha;

alpha = alpha_FW12; % m^2/s
thermal_conductivity = lambda_FW12; %W/m/K
sigma = 5.6695e-8; 
eps = 0.9; 

if dt > 0.5*dx^2/alpha
    disp('bad')
end

%% With flux

%q1 = q_w_stag_TM_new(j); 
n = steps;
T = ones(n,1)*225; %K
b = ones(n,1)*1/2*alpha*dt/(dx^2);
b_mat = 1/2*alpha*dt/(dx^2);

q = [];
for i = 1:steps
    if i == 1
        q(i) = -dx/thermal_conductivity*(qw(i)-sigma*eps*T(1)^4);
    elseif i == steps
        q(i) = -dx/thermal_conductivity*sigma*eps*T(steps)^4;
    else
        q(i) = b_mat*T(i-1)+(1-2*b_mat)*T(i)+b_mat*T(i+1);
    end
end

A = spdiags([-b 1+2*b -b],-1:1,n,n);
A(1,1) = -1;
A(1,2) = 1;
A(steps,steps-1) = -1;
A(steps,steps) = 1;
mat = full(A);

T_front = [];
T_back = [];
T_150 = []; 
T_300 = []; 

for n=1:total_time
    % adjust the right matrix based on the fixed value on boundary
    T = mat\transpose(q);
    T_front(n) = T(1);
    T_150(n) = T(150); 
    T_300(n) = T(300); 
    T_back(n) = T(steps);
    for i = 1:steps
        if i == 1
            q(i) = -dx/thermal_conductivity*(qw(i)-sigma*eps*T(1)^4);
        elseif i == steps
            q(i) = -dx/thermal_conductivity*sigma*eps*T(steps)^4;
        else
            q(i) = b_mat*T(i-1)+(1-2*b_mat)*T(i)+b_mat*T(i+1);
        end
    end
end

 


%%

T_front = [];
T_back = [];

for n=1:total_time
    % adjust the right matrix based on the fixed value on boundary
    T = mat\transpose(q);
    T_front(n) = T(1);
    T_back(n) = T(steps);
    for i = 1:steps
        if i == 1
            q(i) = -dx/thermal_conductivity*(q1-sigma*eps*T(1)^4);
        elseif i == steps
            q(i) = -dx/thermal_conductivity*sigma*eps*T(steps)^4;
        else
            q(i) = b_mat*T(i-1)+(1-2*b_mat)*T(i)+b_mat*T(i+1);
        end
    end
end


%% Plotting

figure(1)
plot(T_front)
hold on
plot(T_back)
legend('front', 'back')

%% Plotting it finally 


