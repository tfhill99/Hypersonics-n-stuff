clear all 
close all

rho_nextel = 2700; % kg/m^3
rho_sigratherm = 92; % kg/m^3
rho_pyrogel = 92; % kg/m^3

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

lambdas_pyrogel = [0.014 0.015 0.020 0.033 0.035 0.40]; 
lambda_pyrogel = mean(lambdas_pyrogel); % W/(K*m)

alpha_nextel = lambda_nextel/(rho_nextel*c_nextel);
alpha_sigratherm = lambda_sigratherm/(rho_sigratherm*c_sigratherm);
alpha_pyrogel = lambda_pyrogel/(rho_pyrogel*c_pyrogel);

