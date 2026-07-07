function J = computeHexJacobianDeterminants(varargin)
% COMPUTEHEXJACOBIANDETERMINANTS Compute center-point Jacobian determinants of a structured hex mesh.
%
% Robust usage:
%   J = computeHexJacobianDeterminants(Domain);
%   J = computeHexJacobianDeterminants(Domain, 'zeroTol', 1e-12);
%   J = computeHexJacobianDeterminants(Xcood_MA, Ycood_MA, Zcood_MA);
%   J = computeHexJacobianDeterminants(Xcood_MA, Ycood_MA, Zcood_MA, 'zeroTol', 1e-12);
%
% Input:
%   Domain or X/Y/Z arrays with size [nely+1, nelx+1, nelz+1].
%
% Name-value options:
%   'autoOrient'  true/false. If true, infer a global orientation convention
%                 from the median raw determinant and multiply all values by
%                 that sign. This avoids reporting all elements as negative
%                 solely because of a reversed node-order convention.
%   'zeroTol'     tolerance for detecting near-zero determinants.
%
% Output fields:
%   detJ           orientation-normalized determinant at each element center.
%   rawDetJ        determinant under the built-in node ordering.
%   center         element center coordinates [nely, nelx, nelz, 3].
%   normalizedDetJ detJ divided by median positive detJ.
%   stats          summary statistics and inversion counts.
%
% Reviewer-response note:
%   This version avoids using nnz as a variable name. In MATLAB, naming a
%   variable nnz shadows the built-in nnz() function and can make expressions
%   such as nnz(mask) be interpreted as logical indexing, causing the error:
%   "逻辑索引包含一个超出数组范围的 true 值".

if nargin < 1
    error('computeHexJacobianDeterminants requires either Domain or X/Y/Z inputs.');
end

%% Parse positional inputs robustly.
firstArg = varargin{1};
if isstruct(firstArg)
    Domain = firstArg;
    requiredFields = {'Xcood_MA','Ycood_MA','Zcood_MA'};
    for k = 1:numel(requiredFields)
        if ~isfield(Domain, requiredFields{k})
            error('Domain must contain field %s.', requiredFields{k});
        end
    end
    X = Domain.Xcood_MA;
    Y = Domain.Ycood_MA;
    Z = Domain.Zcood_MA;
    optArgs = varargin(2:end);
else
    if nargin < 3
        error('Provide either a Domain struct or X, Y, and Z coordinate arrays.');
    end
    X = varargin{1};
    Y = varargin{2};
    Z = varargin{3};
    optArgs = varargin(4:end);
end

p = inputParser;
addParameter(p, 'autoOrient', true, @(v) islogical(v) || isnumeric(v));
addParameter(p, 'zeroTol', 1e-12, @(v) isnumeric(v) && isscalar(v) && v >= 0);
parse(p, optArgs{:});
opts = p.Results;
opts.autoOrient = logical(opts.autoOrient);

%% Coordinate array checks.
X = full(double(X));
Y = full(double(Y));
Z = full(double(Z));

if ~isequal(size(X), size(Y), size(Z))
    error('X, Y, and Z coordinate arrays must have the same size.');
end

[nnyNodes, nnxNodes, nnzNodes] = size(X);
nely = nnyNodes - 1;
nelx = nnxNodes - 1;
nelz = nnzNodes - 1;
if nelx < 1 || nely < 1 || nelz < 1
    error('The structured mesh must contain at least one element in each direction.');
end

%% Center derivatives of an 8-node trilinear hexahedron.
% Node order:
% 1=(-,-,-), 2=(+,-,-), 3=(+,+,-), 4=(-,+,-),
% 5=(-,-,+), 6=(+,-,+), 7=(+,+,+), 8=(-,+,+).
dN = [ ...
    -1, -1, -1;
     1, -1, -1;
     1,  1, -1;
    -1,  1, -1;
    -1, -1,  1;
     1, -1,  1;
     1,  1,  1;
    -1,  1,  1] / 8;

rawDetJ = zeros(nely, nelx, nelz);
center = zeros(nely, nelx, nelz, 3);

for ez = 1:nelz
    for ex = 1:nelx
        for ey = 1:nely
            xn = [X(ey,   ex,   ez); X(ey,   ex+1, ez); X(ey+1, ex+1, ez); X(ey+1, ex,   ez); ...
                  X(ey,   ex,   ez+1); X(ey,   ex+1, ez+1); X(ey+1, ex+1, ez+1); X(ey+1, ex,   ez+1)];
            yn = [Y(ey,   ex,   ez); Y(ey,   ex+1, ez); Y(ey+1, ex+1, ez); Y(ey+1, ex,   ez); ...
                  Y(ey,   ex,   ez+1); Y(ey,   ex+1, ez+1); Y(ey+1, ex+1, ez+1); Y(ey+1, ex,   ez+1)];
            zn = [Z(ey,   ex,   ez); Z(ey,   ex+1, ez); Z(ey+1, ex+1, ez); Z(ey+1, ex,   ez); ...
                  Z(ey,   ex,   ez+1); Z(ey,   ex+1, ez+1); Z(ey+1, ex+1, ez+1); Z(ey+1, ex,   ez+1)];

            xyz = [xn, yn, zn];
            Jmat = xyz.' * dN;  % columns: derivatives w.r.t. xi, eta, zeta
            rawDetJ(ey, ex, ez) = det(Jmat);
            center(ey, ex, ez, :) = mean(xyz, 1);
        end
    end
end

%% Orientation normalization and robust statistics.
rawVec = rawDetJ(:);
finiteRawMask = isfinite(rawVec);
rawFinite = rawVec(finiteRawMask);
nonzeroRaw = rawFinite(abs(rawFinite) > opts.zeroTol);

if isempty(nonzeroRaw)
    orientationSign = 1;
else
    orientationSign = sign(median(nonzeroRaw));
    if orientationSign == 0
        orientationSign = 1;
    end
end

if opts.autoOrient
    detJ = orientationSign * rawDetJ;
else
    detJ = rawDetJ;
end

detVec = detJ(:);
finiteDetMask = isfinite(detVec);
detFinite = detVec(finiteDetMask);
positiveDet = detFinite(detFinite > opts.zeroTol);

if isempty(positiveDet)
    medianPositiveDet = NaN;
else
    medianPositiveDet = median(positiveDet);
end
normalizedDetJ = detJ ./ medianPositiveDet;
normVec = normalizedDetJ(:);
finiteNorm = normVec(isfinite(normVec));

stats = struct();
stats.nelx = nelx;
stats.nely = nely;
stats.nelz = nelz;
stats.numElements = numel(detVec);
stats.orientationSign = orientationSign;
stats.zeroTol = opts.zeroTol;
stats.numNonFinite = sum(~finiteRawMask);
stats.numNegativeRaw = sum(rawVec < -opts.zeroTol);
stats.numPositiveRaw = sum(rawVec > opts.zeroTol);
stats.numNonPositive = sum(finiteDetMask & (detVec <= opts.zeroTol));
stats.numPositive = sum(finiteDetMask & (detVec > opts.zeroTol));

if isempty(detFinite)
    stats.minDetJ = NaN;
    stats.maxDetJ = NaN;
    stats.meanDetJ = NaN;
    stats.medianDetJ = NaN;
    stats.stdDetJ = NaN;
else
    stats.minDetJ = min(detFinite);
    stats.maxDetJ = max(detFinite);
    stats.meanDetJ = mean(detFinite);
    stats.medianDetJ = median(detFinite);
    stats.stdDetJ = std(detFinite);
end

if isempty(finiteNorm)
    stats.minNormalizedDetJ = NaN;
    stats.maxNormalizedDetJ = NaN;
else
    stats.minNormalizedDetJ = min(finiteNorm);
    stats.maxNormalizedDetJ = max(finiteNorm);
end

stats.isLocallyBijective = (stats.numNonPositive == 0) && (stats.numNonFinite == 0);

J = struct();
J.detJ = detJ;
J.rawDetJ = rawDetJ;
J.normalizedDetJ = normalizedDetJ;
J.center = center;
J.stats = stats;
J.description = 'Center-point Jacobian determinants of structured physical macro-elements.';
end
