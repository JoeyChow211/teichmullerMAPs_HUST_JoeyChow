function report = plotJacobianDeterminantReport(DomainOrX, resultDir, varargin)
% PLOTJACOBIANDETERMINANTREPORT Plot reviewer-friendly Jacobian determinant figures.
%
% Main purpose:
%   Quantitatively verify that the 3D physical macro-elements are free from
%   inversion by plotting det(J) at all element centers.
%
% Default behavior:
%   - show figure windows directly;
%   - use the ACTUAL determinant det(J) as the main histogram axis;
%   - by default generate TWO histograms: raw det(J) and normalized det(J)/median(det(J));
%   - do not auto-save figures unless explicitly requested.
%
% Recommended usage for reviewer response:
%   report = plotJacobianDeterminantReport(Domain, [], ...
%       'caseName', 'translation_16x10x20', ...
%       'showFigures', true, ...
%       'saveFiles', false);
%
% Optional:
%   'plotNormalizedHistogram', true  -> also generate det(J)/median(det(J)) histogram.
%   'useLatex', true                 -> use LaTeX interpreters.

p = inputParser;
addParameter(p, 'Y', [], @(v) isnumeric(v) || isempty(v));
addParameter(p, 'Z', [], @(v) isnumeric(v) || isempty(v));
addParameter(p, 'caseName', 'case', @(s) ischar(s) || isstring(s));
addParameter(p, 'showFigures', true, @(v) islogical(v) || isnumeric(v));
addParameter(p, 'numBins', 30, @(v) isnumeric(v) && isscalar(v) && v >= 5);
addParameter(p, 'zeroTol', 1e-12, @(v) isnumeric(v) && isscalar(v) && v >= 0);
addParameter(p, 'makeHistogram', true, @(v) islogical(v) || isnumeric(v));
addParameter(p, 'makeContour', false, @(v) islogical(v) || isnumeric(v));
addParameter(p, 'makeScatter', false, @(v) islogical(v) || isnumeric(v));
addParameter(p, 'plotNormalizedHistogram', true, @(v) islogical(v) || isnumeric(v));
addParameter(p, 'saveFiles', false, @(v) islogical(v) || isnumeric(v));
addParameter(p, 'useLatex', false, @(v) islogical(v) || isnumeric(v));
parse(p, varargin{:});
opts = p.Results;

opts.showFigures = logical(opts.showFigures);
opts.makeHistogram = logical(opts.makeHistogram);
opts.makeContour = logical(opts.makeContour);
opts.makeScatter = logical(opts.makeScatter);
opts.plotNormalizedHistogram = logical(opts.plotNormalizedHistogram);
opts.saveFiles = logical(opts.saveFiles);
opts.useLatex = logical(opts.useLatex);

if nargin < 2 || isempty(resultDir)
    resultDir = fullfile(pwd, 'results', 'reviewer_jacobian');
end
if opts.saveFiles && ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

if isstruct(DomainOrX)
    J = computeHexJacobianDeterminants(DomainOrX, 'zeroTol', opts.zeroTol);
else
    if isempty(opts.Y) || isempty(opts.Z)
        error('When the first input is X, pass Y and Z as name-value arguments.');
    end
    J = computeHexJacobianDeterminants(DomainOrX, opts.Y, opts.Z, 'zeroTol', opts.zeroTol);
end

caseName = char(opts.caseName);
textInterpreter = ternary(opts.useLatex, 'latex', 'tex');
figVis = ternary(opts.showFigures, 'on', 'off');

report = struct();
report.J = J;
report.resultDir = resultDir;
report.histogramFigure = [];
report.contourFigure = [];
report.scatterFigure = [];
report.normalizedHistogramFigure = [];

%% Main histogram using actual determinant values.
if opts.makeHistogram
    detVec = J.detJ(:);
    detVec = detVec(isfinite(detVec));
    fig1 = figure('Visible', figVis, 'Color', 'w', 'Name', 'Jacobian determinant distribution');
    ax1 = axes(fig1); %#ok<LAXES>

    if isempty(detVec)
        text(0.5, 0.5, 'No finite Jacobian values to plot', 'HorizontalAlignment', 'center');
        axis off;
    else
        histogram(ax1, detVec, opts.numBins, ...
            'FaceColor', [0.35, 0.60, 0.80], ...
            'FaceAlpha', 0.90, ...
            'EdgeColor', [1, 1, 1], ...
            'LineWidth', 0.75);
        hold(ax1, 'on');

        minVal = J.stats.minDetJ;
        medVal = J.stats.medianDetJ;
        meanVal = J.stats.meanDetJ;

        xline(ax1, minVal, ':', 'Minimum', ...
            'Color', [0.80, 0.20, 0.20], 'LineWidth', 1.4, ...
            'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'right');
        xline(ax1, medVal, '--', 'Median', ...
            'Color', [0.15, 0.15, 0.15], 'LineWidth', 1.4, ...
            'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
        xline(ax1, meanVal, '-.', 'Mean', ...
            'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.2, ...
            'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');

        xlabel(ax1, 'Jacobian determinant det(J)', 'Interpreter', textInterpreter);
        ylabel(ax1, 'Number of physical macro-elements', 'Interpreter', textInterpreter);
        title(ax1, {sprintf('Jacobian determinant distribution (%d elements)', J.stats.numElements), ...
                    'Element-center values in the physical macro mesh'}, ...
              'Interpreter', textInterpreter, 'FontWeight', 'bold');

        % Make the axis a bit roomier.
        xl = xlim(ax1);
        xr = xl(2) - xl(1);
        if xr > 0
            xlim(ax1, [xl(1) - 0.03*xr, xl(2) + 0.06*xr]);
        end

        if J.stats.numNonPositive == 0
            validityLine = sprintf('All %d determinants are positive', J.stats.numElements);
        else
            validityLine = sprintf('Non-positive determinants: %d', J.stats.numNonPositive);
        end

        statsText = sprintf(['Min = %.4e\n', ...
                             'Median = %.4e\n', ...
                             'Mean = %.4e\n', ...
                             'Std = %.4e\n', ...
                             '%s'], ...
                             J.stats.minDetJ, J.stats.medianDetJ, ...
                             J.stats.meanDetJ, J.stats.stdDetJ, validityLine);
        text(ax1, 0.98, 0.96, statsText, 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'FontSize', 10.5, 'BackgroundColor', [1 1 1], 'Margin', 8, ...
            'EdgeColor', [0.75 0.75 0.75]);

        grid(ax1, 'on');
        ax1.GridAlpha = 0.14;
        ax1.MinorGridAlpha = 0.08;
        box(ax1, 'on');
        applyPaperAxesStyle(ax1, textInterpreter);
        hold(ax1, 'off');
    end

    report.histogramFigure = fig1;
    if opts.saveFiles
        savefig(fig1, fullfile(resultDir, [caseName, '_jacobian_histogram.fig']));
    end
end

%% Optional normalized histogram for supplementary use only.
if opts.plotNormalizedHistogram
    normVec = J.normalizedDetJ(:);
    normVec = normVec(isfinite(normVec));
    figN = figure('Visible', figVis, 'Color', 'w', 'Name', 'Normalized Jacobian determinant distribution');
    axN = axes(figN); %#ok<LAXES>

    if isempty(normVec)
        text(0.5, 0.5, 'No finite Jacobian values to plot', 'HorizontalAlignment', 'center');
        axis off;
    else
        histogram(axN, normVec, opts.numBins, ...
            'FaceColor', [0.55, 0.68, 0.45], ...
            'FaceAlpha', 0.90, ...
            'EdgeColor', [1, 1, 1], ...
            'LineWidth', 0.75);
        hold(axN, 'on');
        xline(axN, 1.0, '--', 'Median-normalized reference = 1', ...
            'Color', [0.15, 0.15, 0.15], 'LineWidth', 1.2, ...
            'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
        xline(axN, J.stats.minNormalizedDetJ, ':', 'Minimum', ...
            'Color', [0.80, 0.20, 0.20], 'LineWidth', 1.4, ...
            'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'right');

        xlabel(axN, 'Normalized determinant det(J) / median(det(J))', 'Interpreter', textInterpreter);
        ylabel(axN, 'Number of physical macro-elements', 'Interpreter', textInterpreter);
        title(axN, {sprintf('Normalized Jacobian determinant distribution (%d elements)', J.stats.numElements), ...
                    'Supplementary plot for scale comparison only'}, ...
              'Interpreter', textInterpreter, 'FontWeight', 'bold');
        grid(axN, 'on'); box(axN, 'on');
        axN.GridAlpha = 0.14;
        applyPaperAxesStyle(axN, textInterpreter);
        hold(axN, 'off');
    end

    report.normalizedHistogramFigure = figN;
    if opts.saveFiles
        savefig(figN, fullfile(resultDir, [caseName, '_jacobian_histogram_normalized.fig']));
    end
end

%% Contour plot of minimum det(J) along y.
if opts.makeContour
    contourField = J.detJ;
    contourField(~isfinite(contourField)) = NaN;
    minAcrossY = nan(size(contourField, 2), size(contourField, 3)); % [nelx, nelz]
    for ex = 1:size(contourField, 2)
        for ez = 1:size(contourField, 3)
            col = contourField(:, ex, ez);
            col = col(isfinite(col));
            if ~isempty(col)
                minAcrossY(ex, ez) = min(col);
            end
        end
    end

    fig2 = figure('Visible', figVis, 'Color', 'w', 'Name', 'Layerwise minimum Jacobian determinant');
    ax2 = axes(fig2); %#ok<LAXES>

    if all(isnan(minAcrossY(:)))
        text(0.5, 0.5, 'No finite Jacobian values to plot', 'HorizontalAlignment', 'center');
        axis off;
    else
        contourf(ax2, 1:size(minAcrossY, 1), 1:size(minAcrossY, 2), minAcrossY.', 18, 'LineStyle', 'none');
        axis(ax2, 'tight');
        xlabel(ax2, 'Element index in x direction', 'Interpreter', textInterpreter);
        ylabel(ax2, 'Element index in h direction', 'Interpreter', textInterpreter);
        title(ax2, {'Contour of minimum Jacobian determinant', ...
                    'Minimum element-center det(J) in each x-h column over the y direction'}, ...
              'Interpreter', textInterpreter, 'FontWeight', 'bold');
        cb = colorbar(ax2);
        ylabel(cb, 'Minimum det(J)', 'Interpreter', textInterpreter);
        if strcmpi(textInterpreter, 'latex')
            cb.TickLabelInterpreter = 'latex';
        else
            cb.TickLabelInterpreter = 'tex';
        end
        colormap(ax2, parula(256));
        applyPaperAxesStyle(ax2, textInterpreter);
    end

    report.contourFigure = fig2;
    if opts.saveFiles
        savefig(fig2, fullfile(resultDir, [caseName, '_jacobian_contour_min_y.fig']));
    end
end

%% Optional 3D scatter plot.
if opts.makeScatter
    centerFlat = reshape(J.center, [], 3);
    detVec = J.detJ(:);
    finiteMask = isfinite(detVec) & all(isfinite(centerFlat), 2);

    fig3 = figure('Visible', figVis, 'Color', 'w', 'Name', 'Jacobian determinant scatter');
    ax3 = axes(fig3); %#ok<LAXES>
    if ~any(finiteMask)
        text(0.5, 0.5, 'No finite Jacobian values to plot', 'HorizontalAlignment', 'center');
        axis off;
    else
        scatter3(ax3, centerFlat(finiteMask,1), centerFlat(finiteMask,2), centerFlat(finiteMask,3), ...
            22, detVec(finiteMask), 'filled');
        axis(ax3, 'equal');
        grid(ax3, 'on'); box(ax3, 'on');
        xlabel(ax3, 'x', 'Interpreter', textInterpreter);
        ylabel(ax3, 'y', 'Interpreter', textInterpreter);
        zlabel(ax3, 'h', 'Interpreter', textInterpreter);
        title(ax3, 'Spatial distribution of det(J) at element centers', ...
            'Interpreter', textInterpreter, 'FontWeight', 'bold');
        cb = colorbar(ax3);
        ylabel(cb, 'det(J)', 'Interpreter', textInterpreter);
        view(ax3, 30, 25);
        applyPaperAxesStyle(ax3, textInterpreter);
    end

    report.scatterFigure = fig3;
    if opts.saveFiles
        savefig(fig3, fullfile(resultDir, [caseName, '_jacobian_scatter.fig']));
    end
end

%% Console summary.
fprintf('\nJacobian determinant summary for %s\n', caseName);
fprintf('  Number of elements        : %d\n', J.stats.numElements);
fprintf('  Minimum det(J)            : %.6e\n', J.stats.minDetJ);
fprintf('  Median det(J)             : %.6e\n', J.stats.medianDetJ);
fprintf('  Mean det(J)               : %.6e\n', J.stats.meanDetJ);
fprintf('  Std det(J)                : %.6e\n', J.stats.stdDetJ);
fprintf('  Min det(J)/median(det(J)) : %.6e\n', J.stats.minNormalizedDetJ);
fprintf('  Non-positive count        : %d\n', J.stats.numNonPositive);
fprintf('  Locally bijective         : %s\n', ternary(J.stats.isLocallyBijective, 'true', 'false'));

end

function applyPaperAxesStyle(ax, textInterpreter)
set(ax, 'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'LineWidth', 1.0, ...
        'TickDir', 'out', ...
        'Layer', 'top');
if strcmpi(textInterpreter, 'latex')
    ax.TickLabelInterpreter = 'latex';
else
    ax.TickLabelInterpreter = 'tex';
end
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
