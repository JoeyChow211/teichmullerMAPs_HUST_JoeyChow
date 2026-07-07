function muInterp = interpolateBeltramiLocalSpline(sampleHeight, muSamples, queryHeight)
% INTERPOLATEBELTRAMILOCALSPLINE Locally interpolate Beltrami data in polar form.
% Author: Joey.Chow
%
% Input:
%   sampleHeight - Sample heights for the sectional tesserae.
%   muSamples    - Beltrami coefficients stored as [numFaces x numSections].
%   queryHeight  - Heights to be interpolated.
%
% Output:
%   muInterp - Interpolated Beltrami coefficients [numFaces x numQueries].
%
% Description:
%   The interpolation is carried out on two decoupled quantities:
%   1) the stretch quotient K = (1 + |mu|) / (1 - |mu|), and
%   2) the unwrapped Beltrami angle theta = angle(mu).
%   Both quantities are evaluated with a local cubic Hermite spline so that
%   each query only depends on neighboring sectional data.

sampleHeight = sampleHeight(:);
queryHeight = queryHeight(:);

if size(muSamples, 2) ~= numel(sampleHeight)
    error('The number of Beltrami samples must match the number of sample heights.');
end

muMagnitude = abs(muSamples);
muMagnitude = min(muMagnitude, 1 - 1e-8);

stretchQuotient = (1 + muMagnitude) ./ max(1 - muMagnitude, 1e-8);
muAngle = unwrap(angle(muSamples), [], 2);

stretchInterp = evaluateLocalCubicHermiteSpline(sampleHeight, stretchQuotient.', queryHeight).';
angleInterp = evaluateLocalCubicHermiteSpline(sampleHeight, muAngle.', queryHeight).';

stretchInterp = max(stretchInterp, 1 + 1e-8);
muMagnitudeInterp = (stretchInterp - 1) ./ (stretchInterp + 1);
muMagnitudeInterp = min(max(muMagnitudeInterp, 0), 1 - 1e-8);

muInterp = muMagnitudeInterp .* exp(1i * angleInterp);

end
