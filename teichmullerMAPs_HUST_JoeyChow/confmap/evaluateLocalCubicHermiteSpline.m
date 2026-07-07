function yq = evaluateLocalCubicHermiteSpline(x, y, xq)
% EVALUATELOCALCUBICHERMITESPLINE Evaluate a local cubic Hermite spline.
% Author: Joey.Chow
%
% Input:
%   x  - Strictly increasing sample coordinates [n x 1].
%   y  - Sample values [n x m].
%   xq - Query coordinates [k x 1].
%
% Output:
%   yq - Interpolated values [k x m].
%
% Description:
%   This routine builds a shape-preserving local spline using a PCHIP-style
%   derivative estimate and evaluates the resulting cubic Hermite segments.
%   The interpolation is local: each segment only depends on nearby samples.

x = x(:);
xq = xq(:);

if size(y, 1) ~= numel(x)
    error('The first dimension of y must equal the number of x samples.');
end

if numel(x) < 2
    error('At least two samples are required for interpolation.');
end

if any(diff(x) <= 0)
    error('Sample coordinates x must be strictly increasing.');
end

if numel(x) == 2
    slope = (y(2, :) - y(1, :)) / (x(2) - x(1));
    yq = y(1, :) + (xq - x(1)) * slope;
    return;
end

h = diff(x);
delta = diff(y, 1, 1) ./ h;
derivatives = computeLocalPchipDerivatives(x, y, h, delta);

yq = zeros(numel(xq), size(y, 2));

for iq = 1:numel(xq)
    if xq(iq) <= x(1)
        idx = 1;
    elseif xq(iq) >= x(end)
        idx = numel(x) - 1;
    else
        idx = find(x <= xq(iq), 1, 'last');
        if idx == numel(x)
            idx = numel(x) - 1;
        end
    end

    hi = x(idx + 1) - x(idx);
    t = (xq(iq) - x(idx)) / hi;

    h00 = 2 * t^3 - 3 * t^2 + 1;
    h10 = t^3 - 2 * t^2 + t;
    h01 = -2 * t^3 + 3 * t^2;
    h11 = t^3 - t^2;

    yq(iq, :) = h00 * y(idx, :) ...
               + h10 * hi * derivatives(idx, :) ...
               + h01 * y(idx + 1, :) ...
               + h11 * hi * derivatives(idx + 1, :);
end

end

function derivatives = computeLocalPchipDerivatives(x, y, h, delta)
% LOCAL_PCHIP_DERIVATIVES Estimate local slopes for the cubic Hermite spline.
% Author: Joey.Chow

n = numel(x);
numSeries = size(y, 2);
derivatives = zeros(n, numSeries);

for j = 1:numSeries
    del = delta(:, j);
    d = zeros(n, 1);

    for i = 2:(n - 1)
        if del(i - 1) == 0 || del(i) == 0 || sign(del(i - 1)) ~= sign(del(i))
            d(i) = 0;
        else
            w1 = 2 * h(i) + h(i - 1);
            w2 = h(i) + 2 * h(i - 1);
            d(i) = (w1 + w2) / (w1 / del(i - 1) + w2 / del(i));
        end
    end

    d(1) = ((2 * h(1) + h(2)) * del(1) - h(1) * del(2)) / (h(1) + h(2));
    if sign(d(1)) ~= sign(del(1))
        d(1) = 0;
    elseif sign(del(1)) ~= sign(del(2)) && abs(d(1)) > abs(3 * del(1))
        d(1) = 3 * del(1);
    end

    d(n) = ((2 * h(end) + h(end - 1)) * del(end) - h(end) * del(end - 1)) / (h(end) + h(end - 1));
    if sign(d(n)) ~= sign(del(end))
        d(n) = 0;
    elseif sign(del(end)) ~= sign(del(end - 1)) && abs(d(n)) > abs(3 * del(end))
        d(n) = 3 * del(end);
    end

    derivatives(:, j) = d;
end

end
