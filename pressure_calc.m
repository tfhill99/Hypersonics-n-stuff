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


thetas = reshape(thetas, [16348, 1]); 
Cp = reshape(Cp, [16348, 1]); 

%% Plotting
f = figure;
trisurf(TR, Cp, 'EdgeColor', 'none');  
colorbar;
