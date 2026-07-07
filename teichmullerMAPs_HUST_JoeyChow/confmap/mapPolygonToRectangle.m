function map = mapPolygonToRectangle(v,f,corner)
% MAPPOLYGONTORECTANGLE Map a disk parameterization to a rectangle.
% Author: Joey.Chow
TR = triangulation(f,v);
B = freeBoundary(TR);
edge = splitBoundaryEdges(B,corner);
temp_edge = edge;
temp_corner = [corner;corner(1)];
for i = 1:4
    for j = 1:4
        if sum(ismember(temp_corner(i:i+1),edge{j})) == 4
            temp_edge{i} = edge{j};
        end
    end
end

vertical = [[edge{4};edge{2}],[zeros(length(edge{4}),1);ones(length(edge{2}),1)]];
horizontal = [[edge{1};edge{3}],[zeros(length(edge{1}),1);ones(length(edge{3}),1)]];

map = reconstructMapFromBeltrami(v,f,zeros(length(f),1),vertical,horizontal);
