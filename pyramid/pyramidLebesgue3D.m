% Computes Lebesgue constants for nodal sets on the pyramid
% input: N - order of basis.
% input: r,s,t - coordinates of points
% input: Nsamples/Nlevels - number of initial sampling points and levels of
%        adaptive hillclimbing. 

function maxleb = pyramidLebesgue3D(N, r, s, t, Nsamples, Nlevels)
if nargin ==0
    N = 10;
    [r s t] = pyramidGLLNodes(N);
    Nsamples = 1000;
    Nlevels = 8;
end

V = pyramidBasisBergot3D(N, r, s, t);

lcoeffs = inv(V);

M = Nsamples;

% randomly sample lebesgue function w/"barycentric" pyramid coordinates
% (linear combinations of vertex positions)
L1 = rand(M,1); L2 = rand(M,1); L3 = rand(M,1); L4 = rand(M,1);
L5 = 1-L1-L2-L3-L4; ids = find(L5>0);

L1 = L1(ids); L2 = L2(ids); L3 = L3(ids);
L4 = L4(ids); L5 = L5(ids);

v1 = [-1 -1 0]; v2 = [1 -1 0]; v3 = [1 1 0]; v4 = [-1 1 0]; v5 = [0 0 1];

h = 1;
for cnt=1:Nlevels
    XYZ = L1*v1 + L2*v2 + L3*v3 + L4*v4 + L5*v5;
    
    f = pyramidLebfn3d(N, lcoeffs, XYZ);
    
    [fsort, ids] = sort(-f);
    ids = ids(1:min(length(ids),M/5));
    
    L1 = L1(ids); L2 = L2(ids); L3 = L3(ids);
    L4 = L4(ids); L5 = L5(ids);
    
    % now perturb again
    h = h/2;
    L1new = L1*ones(1, 10) + h*randn(length(L1),10);
    L2new = L2*ones(1, 10) + h*randn(length(L2),10);
    L3new = L3*ones(1, 10) + h*randn(length(L3),10);
    L4new = L4*ones(1, 10) + h*randn(length(L4),10);
    L5new = 1-L1new-L2new-L3new-L4new;
    
    keep = min(100, length(L1));
    
    L1 = [L1(1:keep);L1new(:)];
    L2 = [L2(1:keep);L2new(:)];
    L3 = [L3(1:keep);L3new(:)];
    L4 = [L4(1:keep);L4new(:)];
    L5 = [L5(1:keep);L5new(:)];
    
    ids = find( (L1+L2+L3+L4+L5)<=1 & L1>=0 & L2>=0 & L3>=0 & L4>=0 & L5>=0);
    L1 = L1(ids); L2 = L2(ids); L3 = L3(ids);
    L4 = L4(ids); L5 = L5(ids);
    
    if 1
        disp(sprintf('N = %i, on level %i, max lebeg const = %3.3e',N,cnt,-full(fsort(1))))
    end
end

XYZ = L1*v1 + L2*v2 + L3*v3 + L4*v4 + L5*v5;
f = pyramidLebfn3d(N, lcoeffs, XYZ);

[~, ids] = sort(-f);
ids = ids(1);
maxleb = f(ids(1));
disp(sprintf('\n'))

function maxleb = PyramidLebesgue3D_new(N,r_in,s_in,t_in)

if nargin ==0
    N = 10;
    [r_in s_in t_in] = pyramidGLLNodes(N);
%     Nsamples = 1000;
%     Nlevels = 8;
end

V = pyramidBasisBergot3D(N,r_in,s_in,t_in);
lcoeffs = inv(V);

Nsample = 4*N;
[r s t] = pyramidEquiNodes(Nsample);

tol = 1e-10;

nlevels = 10;
radius = norm([r(1) s(1) t(1)] - [r(2) s(2) t(2)]);

maxleb = 0;
for level = 1:nlevels    
           
    % VDM matrix    
    V = pyramidBasisBergot3D(N,r,s,t);
    
    XYZ = [r(:) s(:) t(:)];
    leb = pyramidLebfn3d(N,lcoeffs,XYZ);   
    [maxlebVec order] = sort(leb,'descend');
    keep = min(100,length(r)); ids = order(1:keep);
    maxleb = max(maxleb,maxlebVec(1));    
    
    if level < nlevels
        Nsamp = 4;
        rad = linspace(-radius,radius,Nsamp);
        [theta phi] = meshgrid(linspace(0,2*pi,Nsamp));
        
        rn = rad(:)*(cos(theta(:)).*sin(phi(:)))';
        sn = rad(:)*(sin(theta(:)).*sin(phi(:)))';
        tn = rad(:)*cos(phi(:))';
        rn = rn(:); sn = sn(:); tn = tn(:);
        
        rsample = [];   ssample = []; tsample = [];
        for j = 1:length(r(ids))
            rsample = [rsample; r(ids(j))+rn];
            ssample = [ssample; s(ids(j))+sn];
            tsample = [tsample; t(ids(j))+tn];
        end
        r = rsample; s = ssample; t = tsample;
        
        % remove nodes outside boundary of pyramid
        bmask = t < -tol; % bottom boundary
        bmask = bmask | (t > 1+r + tol) | (t > 1-r + tol);
        bmask = bmask | (t > 1+s + tol) | (t > 1-s + tol);
        r(bmask) = []; s(bmask) = []; t(bmask) = [];
        
        radius = radius/2;
    end      
    disp(sprintf('On level %i, maxleb = %d',level,maxleb))
end
disp(sprintf('\n'))




function f = pyramidLebfn3d(N, lcoeffs, XYZ)
X = XYZ(:,1); Y = XYZ(:,2); Z = XYZ(:,3);
[r s t] = pyramidEquiNodes3D(N);

%  vdm  = Vandermonde3D(N, X, Y, Z);
vdm  = pyramidBasisBergot3D(N, X, Y, Z);
fn = vdm*lcoeffs;

f = sum(abs(fn'));