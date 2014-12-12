% Warp and Blend nodes on the surface of the pyramid + interior GLL nodes
% for the pyramid.

function [r s t] = pyramidFaceWBNodes3D(N)

tol = 1e-6;

[r s t] = pyramidDuplexWBNodes3D(N); % use alpha optimized

[rg sg tg] = pyramidGLLNodes3D(N);

x =r ; y = s; z = t;
ids = abs(z-min(z))<tol | abs(1-z-x)<tol | abs(1-z+x)<tol | abs(1-z-y)<tol | abs(1-z+y)<tol; 

x =rg ; y = sg; z = tg;
ids_gll = abs(z-min(z))<tol | abs(1-z-x)<tol | abs(1-z+x)<tol | abs(1-z-y)<tol | abs(1-z+y)<tol; 

r(~ids) = rg(~ids_gll);

s(~ids) = sg(~ids_gll);

t(~ids) = tg(~ids_gll);



