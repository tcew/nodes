% Conical construction of equispaced nodes on the pyramid
% optional arg: t, specifies levels at which to place nodes

function [x y z] = pyramidEquiNodes3D(N,t)

if nargin == 1
    t = linspace(0,1,N+1);
end
x = []; y = []; z= [];
for level = 0:N
    a = (1-t(level+1));
    if level < N        
        r1D = linspace(-a,a,N+1-level);
    else
        r1D = 0;
    end    
    [r s] = meshgrid(r1D);
    x = [x; r(:)];
    y = [y; s(:)];    
    z = [z; t(level+1)*ones(size(r(:)))];
end

