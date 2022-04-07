function [area, areas] = get_triangulation_area(TR)
points = TR.Points; 
connections = TR.ConnectivityList; 
area = 0; 
areas = []; 

for i = 1:length(connections)
    vertices = points(connections(i,:),:);
    p1 = vertices(1,:); 
    p2 = vertices(2,:); 
    p3 = vertices(3,:); 

    curr_area = 0.5 * norm(cross(p2-p1, p3-p1)); 
    areas(i) = curr_area; 
    area = area + curr_area;
end
end
