function s=computeTriangleArea(A,B,C)
% COMPUTETRIANGLEAREA Compute the area of a 2D or 3D triangle.
% Author: Joey.Chow
if length(A)==2  %输入三点是二维平面坐标，变成三维
    AB=[B-A 0];
    BC=[C-B 0];
elseif length(A)==3  %输入三点是三维空间坐标
    AB=B-A;
    BC=C-B;
end
    Z=cross(AB,BC);  %叉乘
    s=1/2*norm(Z);  %取模
end
