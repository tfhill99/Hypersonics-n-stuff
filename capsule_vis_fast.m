figure(1)
TR = stlread('CAD_capsule_3.stl');
trisurf(TR)
hold on 
axis equal
P = incenter(TR);
F = faceNormal(TR);
%quiver3(P(:,1),P(:,2),P(:,3), ...
     %F(:,1),F(:,2),F(:,3),2,'color','black');
hold off




