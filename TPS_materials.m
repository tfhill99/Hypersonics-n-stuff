function[alpha_FW12, alpha_Rescor310M, alpha_Rescor311, alpha_Intek1120, alpha_nextel, alpha_sigratherm, alpha_pyrogel, lambda_FW12, lambda_Rescor310M, lambda_Rescor311, lambda_Intek1120, lambda_nextel, lambda_sigratherm, lambda_pyrogel] = TPS_materials()

max_weight = 8.5; % kg
max_surface_temp = 70; % °C

%% Rigid TPS Materials

rho_FW12 = 2900; % kg/m^3
rho_Rescor310M = 800; % kg/m^3
rho_Rescor311 = 800; % kg/m^3
rho_Intek1120 = 6.4; % kg/m^3

c_FW12 = 1500; % J/K*kg
c_Rescor310M = 1457; % J/K*kg
c_Rescor311 = 1457; % J/K*kg
c_Intek1120 = 1009; % J/K*kg

lambda_FW12 = 3.30; % W/(K*m) Up to 2.2 at 1100 °C
lambda_Rescor310M = 0.18; % W/(K*m) at 1650 °C
lambda_Rescor311 = 0.18; % W/(K*m) at 1430 °C
lambda_Intek1120 = 0.046; % W/(K*m) at any °C

max_temp_FW12 = 1700; % °C
max_temp_Rescor310M = 1650; % °C
max_temp_Rescor311 = 1430; % °C
max_temp_Intek1120 = 300; % °C

alpha_FW12 = lambda_FW12/(rho_FW12*c_FW12);
alpha_Rescor310M = lambda_Rescor310M/(rho_Rescor310M*c_Rescor310M);
alpha_Rescor311 = lambda_Rescor311/(rho_Rescor311*c_Rescor311);
alpha_Intek1120 = lambda_Intek1120/(rho_Intek1120*c_Intek1120);

%% Flexible TPS Materials

rho_nextel = 2700; % kg/m^3
rho_sigratherm = 92; % kg/m^3
rho_pyrogel = 112; % kg/m^3

c_nextel = 1046.7; % J/K*kg
c_sigratherm = 700; % J/K*kg
c_pyrogel = 1150; % J/K*kg

min_thick_nextel = 0.00025; % m
min_thick_sigratherm = 0.00015; % m
min_thick_pyrogel = 0.0050; % m 

max_temp_nextel = 1800; % °C
max_temp_sigratherm = 1000; % °C
max_temp_pyrogel = 650; % °C

lambdas_nextel = [0.110 0.138 0.158 0.168 0.185]; 
lambda_nextel = mean(lambdas_nextel); % W/(K*m)

lambdas_sigratherm = [0.050 0.055 0.055 0.080 0.081 0.11 0.18]; 
lambda_sigratherm = mean(lambdas_sigratherm); % W/(K*m)
%lambda_sigratherm = 0.050;

lambdas_pyrogel = [0.014 0.015 0.020 0.033 0.035 0.40]; 
lambda_pyrogel = mean(lambdas_pyrogel); % W/(K*m)
%lambda_pyrogel = 0.014;

alpha_nextel = lambda_nextel/(rho_nextel*c_nextel);
alpha_sigratherm = lambda_sigratherm/(rho_sigratherm*c_sigratherm);
alpha_pyrogel = lambda_pyrogel/(rho_pyrogel*c_pyrogel);

end

