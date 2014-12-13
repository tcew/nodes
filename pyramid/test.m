% demo driver of various nodal sets on the pyramid. uncomment lines to plot
% different nodal sets and their Lebesgue constant. 

N = 3;

[r s t] = pyramidEquiNodes3D(N); % equidistributed nodes
% [r s t] = pyramidGLLNodes3D(N); % conical/Stroud GLL nodes
% [r s t] = pyramidFaceWBNodes3D(N); % interior GLL, surface Warp-and-Blend nodes
% [r s t] = pyramidWBNodes3D(N); % Interpolatory Warp and Blend nodes
% [r s t] = pyramidDuplexWBNodes3D(N); % Duplex Warp and Blend nodes
% [r s t] = pyramidApproxFekete_greedy(N); % approximate Fekete based on greedy 
% [r s t] = pyramidApproxFekete_QR(N); % approximate Fekete based on QR
if 0 % warning: this is slow, esp for high N. change "if 0" to "if 1" to run.  
    [r s t] = pyramidFekete3D(N); % approximate Fekete based on QR
else
    load('precomputedNodes/feketeNodes')
    rst = fekete{N}; r = rst(:,1); s = rst(:,2); t = rst(:,3);
end

% compute Lebesgue constant of selected nodal set
Nsamples = 10000; Nlevels = 12;
lebval = pyramidLebesgue3D(N,r,s,t,Nsamples,Nlevels);
%[r s t w] = pyramidCubature3D(N); % Stroud quadrature rule, modified from John Burkardt

plot3(r,s,t,'.','markersize',24)

title(sprintf('N = %i, Lebesgue constant = %e',N, lebval))