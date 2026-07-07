function report = generateReviewerJacobianFigure()
% GENERATEREVIEWERJACOBIANFIGURE Reproduce the Jacobian validation requested by the reviewer.
% Author: Joey.Chow
%
% This script uses the same macro discretization as the first numerical
% example in the manuscript: 16 x 10 x 20 = 3200 physical macro-elements.
%
% Default behavior:
%   - display two histograms: the ACTUAL Jacobian determinant det(J) and the
%     normalized determinant det(J)/median(det(J));
%   - do not save figures automatically.
%
% Run in MATLAB:
%   report = generateReviewerJacobianFigure();

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
addpath(genpath(fullfile(rootDir, 'DomainMaterialDef')));
addpath(genpath(fullfile(rootDir, 'modelProcess')));
addpath(genpath(fullfile(rootDir, 'confmap')));
addpath(genpath(fullfile(rootDir, 'geometry')));
addpath(genpath(fullfile(rootDir, 'stlTools')));
addpath(genpath(fullfile(rootDir, 'area')));
addpath(genpath(fullfile(rootDir, 'initialModel')));
addpath(genpath(fullfile(rootDir, 'validation')));

opts = struct();
opts.nelx = 16;
opts.nely = 10;
opts.nelz = 20;
opts.microSize = 31;
opts.numofTesserae = 5;
opts.stlBaseName = 'TranS';
opts.showFigures = true;
opts.saveReviewerFiles = false;
opts.doReconstruction = false;
opts.doReviewerJacobianReport = true;
opts.plotNormalizedReviewerHistogram = true;  % generate the normalized histogram as well.
opts.makeReviewerContour = false;              % by default, generate two histograms only.

[Domain, mapping] = runLoftedSurfaceInterpolation(opts); 

resultDir = fullfile(rootDir, 'results', 'reviewer_jacobian');
report = Domain.jacobianReport;
if opts.saveReviewerFiles
    save(fullfile(resultDir, 'reviewer_jacobian_report.mat'), 'report', 'Domain', 'mapping');
end

fprintf('\nJacobian determinant validation finished.\n');
fprintf('Result directory: %s\n', resultDir);
fprintf('Macro-elements: %d\n', report.J.stats.numElements);
fprintf('min det(J): %.6e\n', report.J.stats.minDetJ);
fprintf('median det(J): %.6e\n', report.J.stats.medianDetJ);
fprintf('mean det(J): %.6e\n', report.J.stats.meanDetJ);
fprintf('min det(J)/median(det(J)): %.6e\n', report.J.stats.minNormalizedDetJ);
fprintf('non-positive det(J) count: %d\n', report.J.stats.numNonPositive);
end
