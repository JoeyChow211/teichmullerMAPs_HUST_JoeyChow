
function plotTriMesh(v, f, arg3)
% PLOTTRIMESH Plot a triangular mesh with optional scalar coloring.
% Author: Joey.Chow

if nargin < 3
    patch('Faces', f, 'Vertices', v, 'FaceColor', [0.6, 1.0, 1.0], 'LineWidth', 0.5);
else
    patch('Faces', f, 'Vertices', v, 'FaceColor', 'flat', ...
        'FaceVertexCData', arg3, 'EdgeColor', 'none', 'LineWidth', 0.5);
    colormap('Copper');
    shading interp;
    set(gcf, 'color', 'w');
end

axis equal tight off;
ax = gca;
ax.Clipping = 'off';

end
