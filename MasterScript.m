C_D = 1.4*ones(1,500)';
C_L = zeros(1,500)';
gamma_0 = -1.4;
alpha = 10;
time = linspace(0, 500, 500)';
C_D_change = {};
C_D_change{1} = C_D;
C_L_change = {};
C_L_change{1} = C_L;
plotting = false;
[M1, V, Z, time, rho, P] = convergence_numsoln(C_D, C_L, time, gamma_0);
[C_D, C_L] = pressure_calc(M1, V, Z, alpha, plotting, rho, P);
C_D_change{2} = C_D;
C_L_change{2} = C_L;

i = 2;
while C_D_change{i} - C_D_change{i-1} > 1
    [M1, V, Z, time, rho, P] = convergence_numsoln(C_D, C_L, time, gamma_0);
    [C_D, C_L] = pressure_calc(M1, V, Z, alpha, plotting, rho, P);
    C_D_change{i+1} = C_D;
    C_L_change{i+1} = C_L;
    i = i + 1;
end


