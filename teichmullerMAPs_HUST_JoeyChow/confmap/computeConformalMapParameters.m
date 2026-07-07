
function [mapping, eigenvalueX, eigenvalueY, eigenvalueZ, F1, F2, F3, cornerIndex] = ...
    computeConformalMapParameters(vertices, faces, nelx, nely, cornerLocation)
% COMPUTECONFORMALMAPPARAMETERS Build a rectangular conformal parameterization and sample grid.
% Author: Joey.Chow

[mapping, cornerIndex, ~] = computeRectangularConformalMap(vertices, faces, cornerLocation);

eigenvalue_xPrarmLinespace = linspace(0, 1, nelx + 1);
eigenvalue_yPrarmlinespace = linspace(0, max(mapping(:, 2)), nely + 1);
[param_EigenvalueX, param_EigenvalueY] = meshgrid(eigenvalue_xPrarmLinespace, eigenvalue_yPrarmlinespace);

F1 = scatteredInterpolant(mapping, vertices(:, 1), 'natural');
F2 = scatteredInterpolant(mapping, vertices(:, 2), 'natural');
F3 = scatteredInterpolant(mapping, vertices(:, 3), 'natural');

eigenvalueX = F1(param_EigenvalueX, param_EigenvalueY);
eigenvalueY = F2(param_EigenvalueX, param_EigenvalueY);
eigenvalueZ = F3(param_EigenvalueX, param_EigenvalueY);

visualizeAngleDistortion(vertices, faces, mapping);

end
