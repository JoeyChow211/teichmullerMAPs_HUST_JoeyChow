
function [Domain, mapping] = processLoftedSurfaces(Domain, mapping, stlBaseName, numofTesserae)
% PROCESSLOFTEDSURFACES Process lofted surfaces from sectional STL files.
% Author: Joey.Chow

[Domain, mapping] = parameterizeLoftedSurfaceSections(Domain, mapping, stlBaseName, numofTesserae);

end
