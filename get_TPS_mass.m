function [mass, masses] = get_TPS_mass(thicknesses, tps_type)

mass = 0; 
masses = []; 

%R-TPS Materials
rho_FW12 = 2900; % kg/m^3
rho_Rescor310M = 800; % kg/m^3
rho_Rescor311 = 800; % kg/m^3
rho_Intek1120 = 6.4; % kg/m^3
%F-TPS Materials
rho_nextel = 2700; % kg/m^3
rho_sigratherm = 92; % kg/m^3
rho_pyrogel = 92; % kg/m^3
R_N = 0.272; %m

if tps_type == 'R'
    desnities = [rho_FW12, rho_Rescor310M, rho_Intek1120, rho_FW12];
    disp('rigid')

    TR = stlread('CAD_capsule_3.stl');
    points = TR.Points(find(TR.Points(:,3)==493),:);
    max_h = max(TR.Points(:,3));
    cone_start = 493;
    h = max_h-cone_start;
    left_lim = max(points(:,2));
    %right_lim = min(points(:,2));
    a = left_lim;
    r = (a^2+h^2)/(2*h);
    V = (((pi*h^2)/3)*(3*r-h))*1e-9;

    V_new = 0;
    h_new = 0;
    r_new = 0;

    for i = 1:length(thicknesses)
        thickness = thicknesses(i);
        h_new = h + thickness;
        r_new = r + thickness;
        V_new = (((pi*h_new^2)/3)*(3*r_new-h_new))*1e-9;
        mass(i) = (V_new-V)*desnities(i);
        V = V_new;
        h = h_new;
        r = r_new;
    end

total_mass = sum(mass);
disp(total_mass);

elseif tps_type == 'F'
    desnities = [rho_nextel, rho_sigratherm, rho_pyrogel, rho_nextel];
    disp('flexible')

    r1 = 
    V = (1/3)*pi*(r1^2 + r1*r2 + r2^2)*h;

else 
    error('TPS type not specified')

end


end



