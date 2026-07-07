
function S = buildFaceToVertexInterpolationMatrix(v, f)
% F2V Build the face-to-vertex interpolation matrix.
% Author: Joey.Chow

TR = triangulation(f, v);
ring = vertexAttachments(TR);
nv = size(v, 1);
nf = size(f, 1);

II = cellfun(@times, ring, num2cell(zeros(nv, 1)), 'UniformOutput', false);
II = cell2mat(cellfun(@plus, II, num2cell((1:nv)'), 'UniformOutput', false)')';
JJ = cell2mat(ring')';
avg = cellfun(@length, ring);

S = sparse(II, JJ, ones(length(JJ), 1), nv, nf);
S = sparse(1:nv, 1:nv, 1 ./ avg) * S;

end
