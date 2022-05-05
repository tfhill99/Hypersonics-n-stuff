function [mass, total_mass] = get_TPS_mass(thicknesses, tps_type)
%thicknesses = [0.005 0.03 0.02 0.005]; %12.9974kg for rigid (what the grads have for optimized...)
%thicknesses = [0.000656 0.010 0.070	0.000656]; %@ time = 500s rigid optimized
%thicknesses = [0.00027 0.006 0.020 0.00025]; % 6.1516 @time = 500s flexible optimized
%thicknesses = [0.00025 0.005 0.0174 0.00025]; % 5.4980 kg for flexible, not converged
%tps_type = 'R';

mass = []; 

%R-TPS Materials
rho_FW12 = 2900; % kg/m^3
rho_Rescor310M = 800; % kg/m^3
rho_Rescor311 = 800; % kg/m^3
rho_Intek1120 = 6.4; % kg/m^3
%F-TPS Materials
rho_nextel = 2700; % kg/m^3
rho_sigratherm = 92; % kg/m^3
rho_pyrogel = 112; % kg/m^3
R_N = 0.272; %m
R_B = 0.765; %m


if tps_type == 'R'
    densities = [rho_FW12, rho_Rescor310M, rho_Intek1120, rho_FW12];
    disp('rigid')

    a = R_N; %left_lim;
    %r = (a^2+h^2)/(2*h);
    h = 0.0692; % measured from CAD
    %V = (((pi*h^2)/3)*(3*r-h))*1e-9;
    V = (1/6)*pi*h*(3*a^2 + h^2);

    V_new = 0;
    h_new = 0;
    a_new = 0;

    flip(thicknesses);
    for i = 1:length(thicknesses)
        thickness = thicknesses(i);
        h_new = h + thickness;
        a_new = a + thickness;
        V_new = (1/6)*pi*h_new*(3*a_new^2 + h_new^2);
        mass(i) = (V_new-V)*densities(i);
        V = V_new;
        h = h_new;
        a = a_new;
    end

total_mass = sum(mass);
disp(total_mass);

elseif tps_type == 'F'
    densities = [rho_nextel, rho_sigratherm, rho_pyrogel, rho_nextel];
    disp('flexible')
    r1 = R_B;
    r2 = R_N;
    h = 0.625; %m
    V = (1/3)*pi*(r1^2 + r1*r2 + r2^2)*h;

    V_new = 0;
    r1_new = 0;
    r2_new = 0;

    flip(thicknesses);
    for i = 1:length(thicknesses)
        thickness = thicknesses(i)*cos(pi/4);
        r1_new = r1 + thickness;
        r2_new = r2 + thickness;
        V_new = (1/3)*pi*(r1_new^2 + r1_new*r2_new + r2_new^2)*h;
        mass(i) = (V_new-V)*densities(i);
        V = V_new;
        r1 = r1_new;
        r2 = r2_new;
    end

total_mass = sum(mass);
disp(total_mass);


else 
    error('TPS type not specified')

end


end
