C_D = 1.4*ones(1,1000);
C_L = zeros(1,1000);
gamma_0 = -1.4;
% need to fix this parameter maybe
time = linspace(0, 500, 1000);
C_D_change = {};
C_D_change{1} = C_D;
C_L_change = {};
C_L_change{1} = C_L;
[M1, V, Z, time] = convergence_numsoln(C_D, C_L, time, gamma_0);
[C_D, C_L] = pressure_calc(M1, V, Z);
C_D_change{2} = C_D;
C_L_change{2} = C_L;

i = 2;
while C_D_change{i} - C_D_change{i-1} > 1
    [M1, V, Z, time] = convergence_numsoln(C_D, C_L, time, gamma_0);
    [C_D, C_L] = pressure_calc(M1, V, Z);
    C_D_change{i+1} = C_D;
    C_L_change{i+1} = C_L;
    i = i + 1;
end


