function A = computeGeneralizedLaplacian(v,f,mu)
% COMPUTEGENERALIZEDLAPLACIAN Compute the generalized Laplacian associated with a Beltrami field.
% Author: Joey.Chow
% Compute the generalized Laplacian.
a_f = (1-2*real(mu)+abs(mu).^2)./(1.0-abs(mu).^2);
b_f = -2*imag(mu)./(1.0-abs(mu).^2);
g_f = (1+2*real(mu)+abs(mu).^2)./(1.0-abs(mu).^2);
f_0 = f(:,1); f_1 = f(:,2); f_2 = f(:,3);
uxv0 = v(f_1,2) - v(f_2,2);
uyv0 = v(f_2,1) - v(f_1,1);
uxv1 = v(f_2,2) - v(f_0,2);
uyv1 = v(f_0,1) - v(f_2,1); 
uxv2 = v(f_0,2) - v(f_1,2);
uyv2 = v(f_1,1) - v(f_0,1);
l = [sqrt(sum(uxv0.^2 + uyv0.^2,2)), sqrt(sum(uxv1.^2 + uyv1.^2,2)), sqrt(sum(uxv2.^2 + uyv2.^2,2))];
s = sum(l,2)*0.5;
area = sqrt(s.*(s-l(:,1)).*(s-l(:,2)).*(s-l(:,3)));
v_00 = (a_f.*uxv0.*uxv0 + 2*b_f.*uxv0.*uyv0 + g_f.*uyv0.*uyv0)./area;
v_11 = (a_f.*uxv1.*uxv1 + 2*b_f.*uxv1.*uyv1 + g_f.*uyv1.*uyv1)./area;
v_22 = (a_f.*uxv2.*uxv2 + 2*b_f.*uxv2.*uyv2 + g_f.*uyv2.*uyv2)./area;
v_01 = (a_f.*uxv1.*uxv0 + b_f.*uxv1.*uyv0 + b_f.*uxv0.*uyv1 + g_f.*uyv1.*uyv0)./area;
v_12 = (a_f.*uxv2.*uxv1 + b_f.*uxv2.*uyv1 + b_f.*uxv1.*uyv2 + g_f.*uyv2.*uyv1)./area;
v_20 = (a_f.*uxv0.*uxv2 + b_f.*uxv0.*uyv2 + b_f.*uxv2.*uyv0 + g_f.*uyv0.*uyv2)./area;

I = [f_0;f_1;f_2;f_0;f_1;f_1;f_2;f_2;f_0];
J = [f_0;f_1;f_2;f_1;f_0;f_2;f_1;f_0;f_2];
V = [v_00;v_11;v_22;v_01;v_01;v_12;v_12;v_20;v_20]/2;
A = sparse(I,J,-V);

end
