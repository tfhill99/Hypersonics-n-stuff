function [C_D, C_L] = pressure_calc(M1, velocity, Z, alpha, plotting, density, pressure, file)
% constants
gamma = 1.4;
P2_P1 = 1 + 2*gamma/(gamma+1)*(M1.^2-1);
Pstag_P2 = ((gamma+1)^2*M1.^2./(4*gamma*M1.^2-2*(gamma-1))).^3;
C_p0 = 2./(gamma*M1.^2).*(P2_P1 .* Pstag_P2-1);
%C_p0_limit =  (4/(gamma+1))*((gamma+1)^2/(4*gamma))^((gamma)/(gamma-1));

% read the mesh and find the vectors
[TR, P, F] = read_mesh(file); 

vx = velocity.*sin(alpha); 
vy = 0; 
vz = velocity.*cos(alpha);  
V = ones(3,length(velocity));
V(1,:) = vx;
V(2,:) = vy;
V(3,:) = vz;
thetas = []; 
Cps = [];

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

thetas = reshape(thetas, [length(thetas), 1]); 
Cps = reshape(Cps, [length(Cps), length(M1)]); 

% finding cps along the velocity vector line
vect_indices_1 = find(P(:,2) > -5); % the y coordinate will be 0 in the x-z plane, where the velocity is applied
vect_indices_2 = find(P(:,2) < 5); 
[vect_indices, eh] = intersect(vect_indices_1, vect_indices_2); 

cross_sec_cp = []; 
cross_sec_X = [];

enter_index = 1; 

for j = 1:length(M1)
    for i = 1:length(vect_indices)
        k = vect_indices(i);
        cross_sec_X(enter_index) = P(k, 1); 
        cross_sec_cp(enter_index, j) = Cps(k, j); 
        enter_index = enter_index + 1; 
    end
end

cross_sec_cp = reshape(cross_sec_cp, [length(cross_sec_cp), length(M1)]);
cross_sec_X = reshape(cross_sec_X, [length(cross_sec_X), 1]);

% write to csv
 D = [cross_sec_X cross_sec_cp]; 
 csvwrite('Cp_along_xzplane', D);

 if plotting == true
    % Cps vs Mach
    Cpps = Cps(1,:);
     
    figure(1)
    axis([0 2 0 12000])
    plot(M1, Cpps,'-.')
    axis([0 12000 0 2])
    legend(' ');
    ylabel('Cp');
    xlabel('Mach');
    title('Cp vs Mach');
 end

 % calculating pressure at each point for each mach
 rho = density(1:length(M1)); 
 P = pressure(1:length(M1)); 
 pressures = [];

 [area, areas] = get_triangulation_area(TR); 

 for i = 1:length(M1)
     pressures(:,i) = Cps(:,i) * 0.5 * rho(i) * velocity(i)^2 + P(i); 
 end

 f_z = []; 
 f_x = []; 

 for i = 1:length(M1)
     f_z_curr = 0; 
     f_x_curr = 0; 
     for j = 1:length(F)
         normals = F(j,:); 
         f_z_curr = f_z_curr + (-1 * pressures(j, i)) * normals(3) * areas(j); 
         f_x_curr = f_x_curr + (-1 * pressures(j, i)) * normals(1) * areas(j);
     end
     f_z(i) = f_z_curr; 
     f_x(i) = f_x_curr; 
 end

 max_R = 765/1000; %m
 reff_area = pi * max_R * max_R; 

 C_L = f_z ./ (rho .* velocity.^2);  
 C_L = C_L / (0.5 * reff_area); 
 C_D = f_x ./ (rho .* velocity.^2);
 C_D = C_D / (0.5 * reff_area);

end


