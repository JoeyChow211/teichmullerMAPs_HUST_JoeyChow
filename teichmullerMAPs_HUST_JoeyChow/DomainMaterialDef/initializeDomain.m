
function Domain = initializeDomain(nelx, nely, nelz, microSize, modelType)
% INITIALIZEDOMAIN Initialize the design-domain data structure.
% Author: Joey.Chow
%
% Input:
%   nelx, nely, nelz - Number of macro cells in each direction.
%   microSize        - Number of nodes per micro-cell direction.
%   modelType        - Integer flag for the geometry type.
%
% Output:
%   Domain - Initialized domain structure.

Domain.nelx = nelx;
Domain.nely = nely;
Domain.nelz = nelz;

Domain.Xcood_MA = zeros(nely + 1, nelx + 1, nelz + 1);
Domain.Ycood_MA = zeros(nely + 1, nelx + 1, nelz + 1);
Domain.Zcood_MA = zeros(nely + 1, nelx + 1, nelz + 1);

Domain.Xcood_MA_MI = zeros(nely * microSize, nelx * microSize, nelz + 1);
Domain.Ycood_MA_MI = zeros(nely * microSize, nelx * microSize, nelz + 1);
Domain.Zcood_MA_MI = zeros(nely * microSize, nelx * microSize, nelz + 1);

Domain.Xcood_MI = zeros(microSize * nely, microSize * nelx, nelz * microSize);
Domain.Ycood_MI = zeros(microSize * nely, microSize * nelx, nelz * microSize);
Domain.Zcood_MI = [];

Domain.modelType = modelType;
Domain.microSize = microSize;

end
