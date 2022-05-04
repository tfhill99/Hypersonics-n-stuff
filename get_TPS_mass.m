function [mass, total_mass] = get_TPS_mass(thicknesses, tps_type)
%thickness = [0.0005 0.0098 0.0459 0.0005]; % Rigid optimized
% produces mass of 1.8519kg w/out safety factor, thickness of 0.0567
%thicknesses = [0.00025 0.0092 0.028 0.00032]; % Flexible optimized
%thicknesses = [0.00035 0.01 0.001 0.00025]; %within mass 5.4550
%thicknesses = [0.0004 0.005 0.004 0.00025];
%produces mass of 11.6752kg w/out safety factor
%thicknesses = [0.003 0.04 0.025 .002]; 11.3254kg (what the grads have LOL)
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

% TR = stlread('CAD_capsule_3.stl');
% points = TR.Points(find(TR.Points(:,3)==493),:);
% max_h = max(TR.Points(:,3));
% cone_start = 493;
% h = max_h-cone_start;
% left_lim = max(points(:,2));
%right_lim = min(points(:,2));

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

    for i = 1:length(thicknesses)
        thickness = thicknesses(i);
        h_new = h + thickness;
        a_new = a + thickness;
        %V_new = (((pi*h_new^2)/3)*(3*r_new-h_new))*1e-9;
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

    for i = 1:length(thicknesses)
        thickness = thicknesses(i);
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
%% calculate mass of rigid
% SA = 4.3565e+06/1e+06 - 1.8385;
% a = sqrt(0.07*(2*R_N-0.07));
% CapSA =  pi*(a^2 + 0.07^2);
% %ConeSA = SA - CapSA;
% rho_rigid = [2900 800 6.4 2900];
% mass_rigid = [];
% %mass_flexible = [];
% for i = 1:length(thickness)
%    mass_rigid(i) = thickness(i)*rho_rigid(i)*CapSA;
%    %mass_flexible(i) = thickness_flexible(i)*rho_flexible(i)*ConeSA;
% end
% mass_sum_rigid = sum(mass_rigid);
