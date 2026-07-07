function mu = computeBeltramiCoefficient(v, f, mapping)
% COMPUTEBELTRAMICOEFFICIENT Compute the Beltrami coefficient of a mapping.
% Author: Joey.Chow
% Compute the Beltrami coefficient of a mapping.
% Input:
% v: nv x (1/2/3) vertex coordinates 
% f: nf x 3 triangulations 
nf = length(f);
M_i = reshape([1:nf;1:nf;1:nf], [1,3*nf]);
M_j = reshape(f', [1,3*nf]);
e_1 = v(f(:,3),1:2) - v(f(:,2),1:2);
e_2 = v(f(:,1),1:2) - v(f(:,3),1:2);
e_3 = v(f(:,2),1:2) - v(f(:,1),1:2);
Area = (-e_2(:,1).*e_1(:,2) + e_1(:,1).*e_2(:,2))'/2;
Area = [Area;Area;Area];

M_x = reshape([e_1(:,2),e_2(:,2),e_3(:,2)]'./Area /2 , [1, 3*nf]);
M_y = -reshape([e_1(:,1),e_2(:,1),e_3(:,1)]'./Area /2 , [1, 3*nf]);
D_x = sparse(M_i,M_j,M_x);
D_y = sparse(M_i,M_j,M_y);

if size(mapping,2) == 2
    Dz = (D_x - 1i*D_y) / 2; Dc = (D_x + 1i*D_y) / 2;
    mu = (Dc*complex(mapping(:,1),mapping(:,2))) ./ (Dz*complex(mapping(:,1),mapping(:,2)));
else
    dX_du = D_x*mapping(:,1);dX_dv = D_y*mapping(:,1);
    dY_du = D_x*mapping(:,2);dY_dv = D_y*mapping(:,2);
    dZ_du = D_x*mapping(:,3);dZ_dv = D_y*mapping(:,3);
    EE = dX_du.^2 + dY_du.^2 + dZ_du.^2;
    GG = dX_dv.^2 + dY_dv.^2 + dZ_dv.^2;
    FF = dX_du.*dX_dv + dY_du.*dY_dv + dZ_du.*dZ_dv;
    mu = (EE - GG + 2 * 1i * FF) ./ (EE + GG + 2*sqrt(EE.*GG - FF.^2));
end
end
