function [TR, P, F] = read_mesh(name)
TR = stlread(name);
P = incenter(TR);
F = faceNormal(TR);
end