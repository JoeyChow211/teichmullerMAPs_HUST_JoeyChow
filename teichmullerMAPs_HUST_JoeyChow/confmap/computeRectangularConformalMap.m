
function [mapping, cornerIndex, TR] = computeRectangularConformalMap(vertices, faces, cornerLocation)
% COMPUTERECTANGULARCONFORMALMAP Map an open simply connected surface to a rectangle.
% Author: Joey.Chow

if size(vertices, 1) < size(vertices, 2)
    vertices = vertices';
end

if size(faces, 1) < size(faces, 2)
    faces = faces';
end

numVertices = size(vertices, 1);
if size(vertices, 2) < 3
    vertices = [vertices, zeros(numVertices, 1)];
end

TR = triangulation(faces, vertices);
triangulationBoundary = freeBoundary(TR);
boundaryPointsIndex = triangulationBoundary(:, 1);

distance = zeros(length(boundaryPointsIndex), 4);
for i = 1:length(boundaryPointsIndex)
    for j = 1:4
        distance(i, j) = sum((TR.Points(boundaryPointsIndex(i), :) - cornerLocation(j, :)).^2);
    end
end

cornerIndexinBoundary = zeros(4, 1);
for i = 1:4
    [~, cornerIndexinBoundary(i)] = min(distance(:, i));
end
cornerIndex = boundaryPointsIndex(cornerIndexinBoundary);

showTriMesh(TR.Points, TR.ConnectivityList);
hold on;
plot3(TR.Points(cornerIndex(1), 1), TR.Points(cornerIndex(1), 2), TR.Points(cornerIndex(1), 3), 'ro', 'MarkerFaceColor', 'r');
plot3(TR.Points(cornerIndex(2), 1), TR.Points(cornerIndex(2), 2), TR.Points(cornerIndex(2), 3), 'ro', 'MarkerFaceColor', 'g');
plot3(TR.Points(cornerIndex(3), 1), TR.Points(cornerIndex(3), 2), TR.Points(cornerIndex(3), 3), 'ro', 'MarkerFaceColor', 'b');
plot3(TR.Points(cornerIndex(4), 1), TR.Points(cornerIndex(4), 2), TR.Points(cornerIndex(4), 3), 'ro', 'MarkerFaceColor', 'y');
title('Initial Surface and Corner Points');
view([-45, 30]);
hold off;

disk = computeDiskConformalMap(TR.Points, TR.ConnectivityList);
rect = mapPolygonToRectangle(disk, TR.ConnectivityList, cornerIndex);

sol_u = rect(:, 1);
sol_v = rect(:, 2);
mu = zeros(size(TR.ConnectivityList, 1), 1);
E_handle = @(h) computeHeightScalingEnergy(h, sol_u, sol_v, TR.Points, TR.ConnectivityList, mu);
h_opt = fminbnd(E_handle, 0, 10);

mapping = [sol_u, h_opt * sol_v];

showTriMesh([mapping(:, 2), mapping(:, 1)], TR.ConnectivityList);
axis equal;

end

function E = computeHeightScalingEnergy(h, sol_u, sol_v, vertices, faces, mu)
% H_ENERGY Auxiliary objective used to balance rectangle height.
% Author: Joey.Chow

h_mu = computeBeltramiCoefficient([sol_u, h * sol_v], faces, vertices);
E = sqrt(sum(abs(h_mu(:, 1) - mu).^2));

end
