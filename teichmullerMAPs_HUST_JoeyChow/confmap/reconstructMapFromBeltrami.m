function map = reconstructMapFromBeltrami(v,f,mu,constraint_x,constraint_y)
% RECONSTRUCTMAPFROMBELTRAMI Reconstruct a quasi-conformal map from boundary constraints.
% Author: Joey.Chow
% quasi-conformal map reconstruction
A = computeGeneralizedLaplacian(v,f,mu);
B = A;
d1 = -A(:,constraint_x(:,1))*constraint_x(:,2);
d1(constraint_x(:,1)) = constraint_x(:,2);
A(constraint_x(:,1),:) = 0;
A(:,constraint_x(:,1)) = 0;
A = A + sparse(constraint_x(:,1),constraint_x(:,1),ones(length(constraint_x(:,1)),1), size(A,1), size(A,2));
umap = A\d1;
d2 = -B(:,constraint_y(:,1))*constraint_y(:,2);
d2(constraint_y(:,1)) = constraint_y(:,2);
B(constraint_y(:,1),:) = 0;
B(:,constraint_y(:,1)) = 0;
B = B + sparse(constraint_y(:,1),constraint_y(:,1),ones(length(constraint_y(:,1)),1), size(B,1), size(B,2));
vmap = B\d2;
map = [umap,vmap];
end
