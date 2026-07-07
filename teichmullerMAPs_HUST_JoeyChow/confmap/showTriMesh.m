
function showTriMesh(v, f, arg3)
% SHOWTRIMESH Display a triangular mesh with optional scalar coloring.
% Author: Joey.Chow

figure;
if nargin < 3
    patch('Faces', f, 'Vertices', v, 'FaceColor', 'w', ...
        'LineWidth', 0.1, 'FaceAlpha', 1, 'EdgeColor', [0.2, 0.4, 0.6]);
    set(gca, 'FontSize', 30, 'FontName', 'Times');
else
    patch('Faces', f, 'Vertices', v, 'FaceColor', 'flat', ...
        'FaceVertexCData', arg3, 'EdgeColor', 'k');
    colorData = load(fullfile(fileparts(mfilename('fullpath')), 'mycolor.mat'));
    set(gca, 'FontSize', 30, 'FontName', 'Times');
    colormap(flipud(colorData.mycolor));
    axis off;
end

axis equal tight;
ax = gca;
ax.Clipping = 'off';

end
