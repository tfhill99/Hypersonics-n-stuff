%% Constants
clc; 
clear; 
alpha = deg2rad(0); 
velocity = 7396.6; %m/s
gamma = 1.4; 
C_p0 = (4/(gamma+1)) * ((gamma+1)^2/(4*gamma))^((gamma)/(gamma-1)); 

%% Read the mesh and find the vectors
TR = stlread('CAD_capsule_3.stl');
P = incenter(TR);
F = faceNormal(TR);
%quiver3(P(:,1),P(:,2),P(:,3), ...
%F(:,1),F(:,2),F(:,3),2,'color','black');

vx = velocity*sin(alpha); 
vy = 0; 
vz = velocity*cos(alpha); 
V = [vx, vy, vz]; 
thetas = []; 
Cp = []; 

for i = 1:length(F)
    N = F(i, :); 
    sin_theta = sum((V.*N)/velocity); 
    thetas(i) = asin(sin_theta); 
    if thetas(i) <= 0
        Cp(i) = 0; 
    else
        Cp(i) = 2 * (sin(thetas(i)))^2;
    end
end


thetas = reshape(thetas, [16348, 1]); 
Cp = reshape(Cp, [16348, 1]); 

%% Read and calculate modified Newtonian
TR = stlread('CAD_capsule_3.stl');
P = incenter(TR);
F = faceNormal(TR);
%quiver3(P(:,1),P(:,2),P(:,3), ...
%F(:,1),F(:,2),F(:,3),2,'color','black');

vx = velocity*sin(alpha); 
vy = 0; 
vz = velocity*cos(alpha); 
V = [vx, vy, vz]; 
thetas = []; 
Cp = []; 

for i = 1:length(F)
    N = F(i, :); 
    sin_theta = sum((V.*N)/velocity); 
    thetas(i) = asin(sin_theta); 
    if thetas(i) <= 0
        Cp(i) = 0; 
    else
        Cp(i) = C_p0 * (sin(thetas(i)))^2;
    end
end


thetas = reshape(thetas, [length(thetas), 1]); 
Cp = reshape(Cp, [length(Cp), 1]); 


%% Finding CP along the velocity application line
vect_indices_1 = find(P(:,2) > -5); % the y coordinate will be 0 in the x-z plane, where the velocity is applied
vect_indices_2 = find(P(:,2) < 5); 
[vect_indices, eh] = intersect(vect_indices_1, vect_indices_2); 

cross_sec_cp = []; 
cross_sec_X = [];

enter_index = 1; 

for i = 1:length(vect_indices)
    k = vect_indices(i);
    cross_sec_X(enter_index) = P(k, 1); 
    cross_sec_cp(enter_index) = Cp(k); 
    enter_index = enter_index + 1; 
end

cross_sec_cp = reshape(cross_sec_cp, [length(cross_sec_cp), 1]);
cross_sec_X = reshape(cross_sec_X, [length(cross_sec_X), 1]);

%% Plotting
figure(1); 
trisurf(TR, Cp, 'EdgeColor', 'none');  
hold on; 
axis equal; 
colorbar;
colormapeditor; 

%% Plotting Cps along cross section
figure(2); 
scatter(cross_sec_X, cross_sec_cp); 
