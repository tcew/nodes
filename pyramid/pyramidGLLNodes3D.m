% Conical construction of GLL nodes on the pyramid

function [x y z] = pyramidGLLNodes3D(N)

t = JacobiGL(0,0,N);
t = (t+1)/2;
x = []; y = []; z= [];
for level = 0:N
    a = 1-t(level+1);
    if level < N
        r1D = a*JacobiGL(0,0,N-level);
%         r1D = a*linspace(-1,1-tol,N+1-level);
    else
        r1D = 0;
    end    
    [r s] = meshgrid(r1D);
    x = [x; r(:)];
    y = [y; s(:)];    
    z = [z; t(level+1)*ones(size(r(:)))];
end
