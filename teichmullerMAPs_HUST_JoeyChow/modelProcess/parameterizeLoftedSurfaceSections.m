
function [Domain, mapping] = parameterizeLoftedSurfaceSections(Domain, mapping, stlBaseName, numofTesserae)
% PARAMETERIZELOFTEDSURFACESECTIONS Parameterize each sectional lofted surface independently.
% Author: Joey.Chow
%
% Input:
%   Domain        - Design-domain structure.
%   mapping       - Cell array for mapping coordinates.
%   stlBaseName   - Common STL filename prefix.
%   numofTesserae - Number of lofted sections.
%
% Output:
%   Domain  - Updated domain data.
%   mapping - Parameterization result for each section.

Xcood_MAinTesserae = zeros(Domain.nely + 1, Domain.nelx + 1, numofTesserae);
Ycood_MAinTesserae = zeros(Domain.nely + 1, Domain.nelx + 1, numofTesserae);
Zcood_MAinTesserae = zeros(Domain.nely + 1, Domain.nelx + 1, numofTesserae);

Domain.Xcood_MA_MIinTesserae = zeros(Domain.nely * Domain.microSize, Domain.nelx * Domain.microSize, numofTesserae);
Domain.Ycood_MA_MIinTesserae = zeros(Domain.nely * Domain.microSize, Domain.nelx * Domain.microSize, numofTesserae);
Domain.Zcood_MA_MIinTesserae = zeros(Domain.nely * Domain.microSize, Domain.nelx * Domain.microSize, numofTesserae);

h_inESpace = [];
Domain.H = zeros(numofTesserae, 1);
Domain.Beltrami_coeff = cell(numofTesserae, 1);
Domain.faceHeight = zeros(numofTesserae, 1);
Domain.vertices = cell(numofTesserae, 1);
Domain.faces = cell(numofTesserae, 1);

for i = 1:numofTesserae
    stlFilename = sprintf('%s%d.stl', stlBaseName, i);
    [vertices, faces, ~, name] = readStlMesh(stlFilename);

    Domain.vertices{i} = vertices;
    Domain.faces{i} = faces;

    h_inESpace = [h_inESpace; vertices(1, 3)];
    plotStlMesh(vertices, faces, name);

    cornerLocation = [8000, -4000, 0;
                      8000,  4000, 0;
                     -8000,  4000, 0;
                     -8000, -4000, 0];

    [mapping{i}, eigenvalueX, eigenvalueY, eigenvalueZ, F1, F2, F3, ~] = ...
        computeConformalMapParameters(vertices, faces, Domain.nelx, Domain.nely, cornerLocation);

    Xcood_MAinTesserae(:, :, i) = eigenvalueX;
    Ycood_MAinTesserae(:, :, i) = eigenvalueY;
    Zcood_MAinTesserae(:, :, i) = eigenvalueZ;

    X_linespace = linspace(0, 1, Domain.microSize * Domain.nelx);
    Y_linespace = linspace(0, max(mapping{i}(:, 2)), Domain.microSize * Domain.nely);
    [x_infill, y_infill] = meshgrid(X_linespace, Y_linespace);

    Domain.Xcood_MA_MIinTesserae(:, :, i) = F1(x_infill, y_infill);
    Domain.Ycood_MA_MIinTesserae(:, :, i) = F2(x_infill, y_infill);
    Domain.Zcood_MA_MIinTesserae(:, :, i) = F3(x_infill, y_infill);

    Domain.H(i) = max(mapping{i}(:, 2));
    Domain.Beltrami_coeff{i} = computeBeltramiCoefficient(vertices, faces, mapping{i});
    Domain.faceHeight(i) = vertices(1, 3);
end

H_Tresserae = linspace(0, 1, numofTesserae);
Z_linspace_param = linspace(0, 1, Domain.nelz + 1);

for i = 1:(Domain.nely + 1)
    for j = 1:(Domain.nelx + 1)
        Domain.Xcood_MA(i, j, :) = makima(H_Tresserae, reshape(Xcood_MAinTesserae(i, j, :), 1, []), Z_linspace_param);
        Domain.Ycood_MA(i, j, :) = makima(H_Tresserae, reshape(Ycood_MAinTesserae(i, j, :), 1, []), Z_linspace_param);
        Domain.Zcood_MA(i, j, :) = makima(H_Tresserae, reshape(Zcood_MAinTesserae(i, j, :), 1, []), Z_linspace_param);
    end
end

for i = 1:(Domain.nely * Domain.microSize)
    for j = 1:(Domain.nelx * Domain.microSize)
        Domain.Xcood_MA_MI(i, j, :) = makima(H_Tresserae, reshape(Domain.Xcood_MA_MIinTesserae(i, j, :), 1, []), Z_linspace_param);
        Domain.Ycood_MA_MI(i, j, :) = makima(H_Tresserae, reshape(Domain.Ycood_MA_MIinTesserae(i, j, :), 1, []), Z_linspace_param);
        Domain.Zcood_MA_MI(i, j, :) = makima(H_Tresserae, reshape(Domain.Zcood_MA_MIinTesserae(i, j, :), 1, []), Z_linspace_param);
    end
end

end
