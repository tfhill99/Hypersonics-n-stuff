function [C_D, C_L, C_M, pressures, areas, f_z, f_x, Cps, C_p0] = pressure_calc(M1, velocity, alpha_init, plotting, density, pressure, file)
% constants
gamma = 1.4;

alpha = deg2rad(alpha_init); 

% nose pressure, contingent on shock theory and M1 dependance
P2_P1 = 1 + 2*gamma/(gamma+1)*(M1.^2-1);
Pstag_P2 = ((gamma+1)^2*M1.^2./(4*gamma*M1.^2-2*(gamma-1))).^3.5;
C_p0 = 2./(gamma*M1.^2).*(P2_P1 .* Pstag_P2-1);

% read the mesh and find the vectors
[TR, Centers, F] = read_mesh(file); 

csvwrite('Triangulation Points', TR.Points); 
csvwrite('Triangulation Connections', TR.ConnectivityList); 

disp('loaded mesh')

vx = velocity.*sin(alpha); 
vy = 0; 
vz = velocity.*cos(alpha);  
V = ones(3,length(velocity));
V(1,:) = vx;
V(2,:) = vy;
V(3,:) = vz;
thetas = []; 
Cps = [];
Cps_nosecone = zeros(length(F), length(M1));
Cps_base = zeros(length(F), length(M1));

disp('starting cp calc')

for j = 1:length(M1)
    for i = 1:length(F)
        N = F(i, :); 
        sin_theta = sum((V(:,j).'.*N)/velocity(j)); 
        thetas(i) = asin(sin_theta); % make 3 rows
        if thetas(i) <= 0
            Cps(i,j) = 0; % Base pressure approximation
            if N(3) < 0
                Cps(i,j) = -2/(gamma * M1(j)^2);
                Cps_base(i, j) = -2/(gamma * M1(j)^2);
            end
            
        else
            Cps(i,j) = C_p0(j) * (sin(thetas(i)))^2; % Newtonian Theory
            Cps_nosecone(i,j) = C_p0(j) * (sin(thetas(i)))^2;
        end
    end
end

disp('ending cp calc')

%thetas = reshape(thetas, [length(thetas), 1]); 
%Cps = reshape(Cps, [length(Cps), length(M1)]); 

 disp('starting pressure calc')
 % already been interpolated in convergence_numsoln.m
 rho = density;
 P = pressure;
 pressures = [];
 pressures_nosecone = [];
 pressures_base = [];

 [~, areas] = get_triangulation_area(TR); 
 
 areas = areas / 1e6; 

 for i = 1:length(M1)
     % Find the pressure for a given Cp (based off Cp definition)
     pressures(:,i) = Cps(:,i) * 0.5 * rho(i) * velocity(i)^2 + P(i);
     pressures_nosecone(:,i) = Cps_nosecone(:,i) * 0.5 * rho(i) * velocity(i)^2 + P(i);
     pressures_base(:,i) = Cps_base(:,i) * 0.5 * rho(i) * velocity(i)^2 + P(i);
 end

 % Force vectors
 f_z = []; 
 f_x = []; 
 f_z_nosecone = [];
 f_x_nosecone = [];
 f_z_base = [];
 f_x_base = [];
 Mys = [];
 Mys_nosecone = [];
 Mys_base = [];
 Centers_X = Centers(:,1);
 Centers_Z = Centers(:,3);
 Centers_X = Centers_X/1000;
 Centers_Z = Centers_Z/1000;
 cg_x = 0;
 cg_z = 0.2809269856; % measurement from CAD

 for i = 1:length(M1)
     f_z_curr = 0; 
     f_x_curr = 0;
     My_curr = 0;
     f_z_curr_nosecone = 0; 
     f_x_curr_nosecone = 0;
     My_curr_nosecone = 0;
     f_z_curr_base = 0; 
     f_x_curr_base = 0;
     My_curr_base = 0;
     
     for j = 1:length(F)
         % find the normal in z and x direction i.e. index 3 and 1
         normals = F(j,:); 
         My_curr = My_curr - (- (Centers_X(j)-cg_x)* (-pressures(j, i) * normals(3)) + (Centers_Z(j)-cg_z)* (-pressures(j, i))* normals(1)) * areas(j);
         My_curr_nosecone = My_curr_nosecone - (- (Centers_X(j)-cg_x)* (-pressures_nosecone(j, i) * normals(3)) + (Centers_Z(j)-cg_z)* (-pressures_nosecone(j, i)) * normals(1)) * areas(j);
         My_curr_base = My_curr_base - (- (Centers_X(j)-cg_x)* (-pressures_base(j, i) * normals(3)) + (Centers_Z(j)-cg_z)* (-pressures_base(j, i)) * normals(1)) * areas(j);
         f_z_curr = f_z_curr - (pressures(j, i)) * normals(3) * areas(j); 
         f_x_curr = f_x_curr - (pressures(j, i)) * normals(1) * areas(j);
         f_z_curr_nosecone = f_z_curr_nosecone - (pressures_nosecone(j, i)) * normals(3) * areas(j); 
         f_x_curr_nosecone = f_x_curr_nosecone - (pressures_nosecone(j, i)) * normals(1) * areas(j);
         % error in mesh with Normal = 0, causes CL = CL_noscone
         f_z_curr_base = f_z_curr_base - (pressures_base(j, i)) * normals(3) * areas(j); 
         f_x_curr_base = f_x_curr_base - (pressures_base(j, i)) * normals(1) * areas(j);
         %error fix with below line
         %f_x_curr_nosecone = f_x_curr - f_x_curr_base;
         
     end
     f_z(i) = f_z_curr; 
     f_x(i) = f_x_curr;
     f_z_nosecone(i) = f_z_curr_nosecone; 
     f_x_nosecone(i) = f_x_curr_nosecone; 
     f_z_base(i) = f_z_curr_base; 
     f_x_base(i) = f_x_curr_base; 
     Mys(i) = My_curr;
     Mys_nosecone(i) = My_curr_nosecone;
     Mys_base(i) = My_curr_base;
 end

 disp('ending other calc')

D = f_z * cos(alpha) - f_x * sin(alpha); 
L = f_z * sin(alpha) + f_x * cos(alpha);
D_nosecone = f_z_nosecone * cos(alpha) - f_x_nosecone * sin(alpha); 
L_nosecone = f_z_nosecone * sin(alpha) + f_x_nosecone * cos(alpha);  
D_base = f_z_base * cos(alpha) - f_x_base * sin(alpha); 
L_base = f_z_base * sin(alpha) + f_x_base * cos(alpha);

D = -1 * D;
D_nosecone = -1 * D_nosecone;
D_base = -1 * D_base;

 max_R = 765/1000; %m
 reff_area = pi * max_R * max_R; 

 lref = 0.43 * 2 * max_R;

 rho = density;

 C_D = transpose(D) ./ (rho .* velocity.^2);  
 C_D = C_D / (0.5 * reff_area); 
 C_D_nosecone = transpose(D_nosecone) ./ (rho .* velocity.^2);  
 C_D_nosecone = C_D_nosecone / (0.5 * reff_area); 
 C_D_base = transpose(D_base) ./ (rho .* velocity.^2);  
 C_D_base = C_D_base / (0.5 * reff_area); 
 C_L = transpose(L) ./ (rho .* velocity.^2);
 C_L = C_L / (0.5 * reff_area);
 C_L_nosecone = transpose(L_nosecone) ./ (rho .* velocity.^2);
 C_L_nosecone = C_L_nosecone / (0.5 * reff_area);
 C_L_base = transpose(L_base) ./ (rho .* velocity.^2);
 C_L_base = C_L_base / (0.5 * reff_area);
 C_M = transpose(Mys)./(rho .* velocity.^2);
 C_M = C_M / (0.5 * lref * reff_area);
 C_M_nosecone = transpose(Mys_nosecone)./(rho .* velocity.^2);
 C_M_nosecone = C_M_nosecone/ (0.5 * lref * reff_area);
 C_M_base= transpose(Mys_base)./(rho .* velocity.^2);
 C_M_base = C_M_base / (0.5 * lref * reff_area);

%D = [C_L C_L_nosecone C_L_base C_D C_D_nosecone C_D_base C_M C_M_nosecone C_M_base];
%name = "CL_CLnosecone_CLbase_CD_CDnosecone_CDbase_CM_CMnosecone_CMbase_at_" + 0; 
%csvwrite(name, D); 

end
