% Jane Yang 2018
% Modified from https://www.mathworks.com/matlabcentral/answers/93623-how-do-i-plot-the-line-of-intersection-between-two-surfaces
% Goal: Finds the intersection between two surfaces 
% NOTE: Called in visualize_manifolds.m

function [xL, yL, zL]= find_inter(surf1, surf2)
[x4,y4] = meshgrid(linspace(0,4,1000),linspace(0,4,1000));
op1 = surf1(x4,y4);
op2 = surf2(x4,y4);

% Find the difference field.
zdiff = op1 - op2;
% Find the contour where the difference (on the surface) is zero.
[M,c] = contour(x4,y4,zdiff,[0 0]);
% Extract the x- and y-locations from the contour matrix C.
xL = M(1, 2:end);
yL = M(2, 2:end);
% Interpolate on the first surface to find z-locations for the intersection
% line.
zL = interp2(x4, y4, op1, xL, yL);

delete(c)
end