function [C_D, C_L, pressures, areas, f_z, f_x] = pressure_calc(M1, velocity, alpha, plotting, density, pressure, file)
% constants
gamma = 1.4;

% nose pressure, contingent on shock theory and M1 dependance
P2_P1 = 1 + 2*gamma/(gamma+1)*(M1.^2-1);
Pstag_P2 = ((gamma+1)^2*M1.^2./(4*gamma*M1.^2-2*(gamma-1))).^3;
C_p0 = 2./(gamma*M1.^2).*(P2_P1 .* Pstag_P2-1);

% read the mesh and find the vectors
[TR, ~, F] = read_mesh(file); 

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

disp('starting cp calc')

for j = 1:length(M1)
    for i = 1:length(F)
        N = F(i, :); 
        sin_theta = sum((V(:,j).'.*N)/velocity(j)); 
        thetas(i) = asin(sin_theta); % make 3 rows
        if thetas(i) <= 0
            Cps(i,j) = -2/(gamma * M1(j)^2); % Base pressure approximation
        else
            Cps(i,j) = C_p0(j) * (sin(thetas(i)))^2; % Newtonian Theory
        end
    end
end

disp('ending cp calc')

thetas = reshape(thetas, [length(thetas), 1]); 
Cps = reshape(Cps, [length(Cps), length(M1)]); 

 disp('starting pressure calc')
 % already been interpolated in convergence_numsoln.m
 rho = density;
 P = pressure;
 pressures = [];

 [~, areas] = get_triangulation_area(TR); 

 areas = areas / 1e6; 

 for i = 1:length(M1)
     % Find the pressure for a given Cp (based off Cp definition)
     pressures(:,i) = Cps(:,i) * 0.5 * rho(i) * velocity(i)^2 + P(i); 
 end

 % Force vectors
 f_z = []; 
 f_x = []; 

 for i = 1:length(M1)
     f_z_curr = 0; 
     f_x_curr = 0; 
     for j = 1:length(F)
         normals = F(j,:); 
         % find the normal in z and x direction i.e. index 3 and 1
         f_z_curr = f_z_curr - (pressures(j, i)) * normals(3) * areas(j); 
         f_x_curr = f_x_curr - (pressures(j, i)) * normals(1) * areas(j);
     end
     f_z(i) = f_z_curr; 
     f_x(i) = f_x_curr; 
 end

 disp('ending other calc')

 max_R = 765/1000; %m
 reff_area = pi * max_R * max_R; 

 C_D = transpose(f_z) ./ (rho .* velocity.^2);  
 C_D = C_D / (0.5 * reff_area); 
 C_L = transpose(f_x) ./ (rho .* velocity.^2);
 C_L = C_L / (0.5 * reff_area);

end


