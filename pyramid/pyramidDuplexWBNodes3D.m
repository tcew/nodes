% glues nodes from two tets together to get Duplex pyramid nodes
% 
% alphashare* = blending parameters for shared face (edge1-edge3, then face
% blending parameter). If not specified, defaults to alpha.

function [r s t] = pyramidDuplexWBNodes3D(N,alpha,alphaShare1,alphaShare2,alphaShare3,alphaShare4)

tol = 1e-10;

alphastore = [0;0;0;0.1002; 1.1332;1.5608;1.3413;1.2577;1.1603;...
    1.10153;0.6080;0.4523;0.8856;0.8717;0.9655];
if nargin==0
    N = 5;
    
    alpha = alphastore(N);
    alphaShare1 = alphastore(N);
    alphaShare2 = alphastore(N);
    alphaShare3 = alphastore(N);
    alphaShare4 = alphastore(N);
end

% alpha
if nargin == 1
    alpha = alphastore(N);
    alphaShare1 = alpha;
    alphaShare2 = alpha;
    alphaShare3 = alpha;
    alphaShare4 = alpha;
end

% if alpha specified
if nargin ==2
    alphaShare1 = alpha;
    alphaShare2 = alpha;
    alphaShare3 = alpha;
    alphaShare4 = alpha;
end
    
[r s t] = splitTetNodes3D(N,alpha,alphaShare1,alphaShare2,alphaShare3,alphaShare4);

[L1 L2 L3 L4] = rsttobary(r,s,t);

% pyramid vertices
v1 = [-1 -1 0]; v2 = [ 1 -1 0]; v3 = [ 1  1 0]; v4 = [-1  1 0]; v5 = [ 0  0 1];

% first half of pyramid
vv1 = v3; vv2 = v2; vv3 = v4; vv4 = v5;
XYZ = L1*vv1 + L2*vv2 + L3*vv3 + L4*vv4;
r = XYZ(:,1);s = XYZ(:,2);t = XYZ(:,3);

% remove interface nodes
inds = abs(r + s) > tol;
r = r(inds);t = t(inds);s = s(inds);

% second half of pyramid
vv1 = v1; vv2 = v2; vv3 = v4; vv4 = v5;
XYZ = L1*vv1 + L2*vv2 + L3*vv3 + L4*vv4;
r = [r; XYZ(:,1)];s = [s; XYZ(:,2)];t = [t; XYZ(:,3)];

plotflag = nargin==0;
if plotflag
    
    %inds = abs(t)>tol & abs(1-s-t)>tol & abs(1+s-t)>tol & abs(1-r-t)>tol & abs(1+r-t)> tol;
    %r = r(inds);t = t(inds);s = s(inds);
    %     [re se te] = pyramidEquiNodes(N);
    %     plot3(re,se,te,'cs')
    [re se te] = pyramidEquiNodes(1);
    hold on
    inds = [1 2 4 3 1 5 2 5 3 5 4 5];
    re = re(inds);
    se = se(inds);
    te = te(inds);
    for e = 1:length(inds)-1
        plot3(re(e:e+1),se(e:e+1),te(e:e+1),'k-','linewidth',4)
        %     text(r(e),s(e),t(e),num2str(e))
    end    
    
    
    plot3(r,s,t,'k.','markersize',24);
    view(135,0)
    %     axis square
    %     keyboard
end



function [r s t] = splitTetNodes3D(N,alpha,alphaShare1,alphaShare2,alphaShare3,alphaShare4)

if nargin ==0
    N = 3;
    alpha = 0;% choose optimized blending parameter
end

% total number of nodes and tolerance
Np = (N+1)*(N+2)*(N+3)/6; tol = 1e-10;

[r,s,t] = EquiNodes3D(N); % create equidistributed nodes
L1 = (1+t)/2; L2 = (1+s)/2; L3 = -(1+r+s+t)/2; L4 = (1+r)/2;

% set vertices of tetrahedron
v1 = [-1, -1/sqrt(3), -1/sqrt(6)]; v2 = [ 1, -1/sqrt(3),-1/sqrt(6)];
v3 = [ 0,  2/sqrt(3), -1/sqrt(6)]; v4 = [ 0,  0,         3/sqrt(6)];

% orthogonal axis tangents on faces 1-4
t1(1,:) = v2-v1;          t1(2,:) = v2-v1;
t1(3,:) = v3-v2;          t1(4,:) = v3-v1;
t2(1,:) = v3-0.5*(v1+v2); t2(2,:) = v4-0.5*(v1+v2);
t2(3,:) = v4-0.5*(v2+v3); t2(4,:) = v4-0.5*(v1+v3);

for n=1:4 % normalize tangents
    t1(n,:) = t1(n,:)/norm(t1(n,:)); t2(n,:) = t2(n,:)/norm(t2(n,:));
end

% Warp and blend for each face (accumulated in shiftXYZ)
XYZ = L3*v1+L4*v2+L2*v3+L1*v4; % form undeformed coordinates
shift = zeros(size(XYZ));
for face=1:4
    if(face==1); La = L1; Lb = L2; Lc = L3; Ld = L4; end;
    if(face==2); La = L2; Lb = L1; Lc = L3; Ld = L4; end;
    if(face==3); La = L3; Lb = L1; Lc = L4; Ld = L2; end;
    if(face==4); La = L4; Lb = L1; Lc = L3; Ld = L2; end;
    
    blend = Lb.*Lc.*Ld;   % compute volume blending
    
    % compute warp tangential to face
    if face==1
        x = XYZ(:,1); y = XYZ(:,2); z = XYZ(:,3);
        %fids = abs(z-min(z))<tol;
        
        [r s t] = xyztorst(x,y,z);
        warp1r = Warpfactor1D(N,r); warp2s = Warpfactor1D(N,s);
        
        % get warps in terms of xy coordinates
        warp1 = 1*warp1r + warp2s*cos(pi/3);
        warp2 = 0*warp1r + warp2s*sin(pi/3);
        denom = (Lb+.5*La).*(Lc+.5*La).*(Ld+.5*La);   % modify linear blend
        
        ids = find(denom>tol);
        blend(ids) = (1+(alpha.*La(ids)).^2).*blend(ids)./denom(ids);
        
    elseif face==3 % shared face of pyramid
        [warp1 warp2] = WarpShiftFace3D(N, alphaShare1, alphaShare2, alphaShare3, La, Lb, Lc, Ld);
        %         warp1 = 0*warp1; warp2 = 0*warp2;
        %blend = blend.*(1 + (alphaShare*La).^2);
        denom = (Lb+.5*La).*(Lc+.5*La).*(Ld+.5*La);   % modify linear blend
        ids = find(denom>tol);
        blend(ids) = (1+(alphaShare4.*La(ids)).^2).*blend(ids)./denom(ids);
        %        blend(ids) = (1+(alpha.*La(ids)).^2).*blend(ids)./denom(ids);
        
    else
        [warp1 warp2] = WarpShiftFace3D(N, alpha, alpha, alpha, La, Lb, Lc, Ld);
        denom = (Lb+.5*La).*(Lc+.5*La).*(Ld+.5*La);   % modify linear blend
        ids = find(denom>tol);
        blend(ids) = (1+(alpha.*La(ids)).^2).*blend(ids)./denom(ids);
    end
    % compute warp & blend
    shift = shift+(blend.*warp1)*t1(face,:) + (blend.*warp2)*t2(face,:);
    
    % fix face warp
    ids = find(La<tol & ( (Lb>tol) + (Lc>tol) + (Ld>tol) < 3));
    shift(ids,:) = warp1(ids)*t1(face,:) + warp2(ids)*t2(face,:);
    
end
XYZ = XYZ + shift;
X = XYZ(:,1); Y = XYZ(:,2); Z = XYZ(:,3);
[r,s,t] = xyztorst(X,Y,Z);
% leb = Lebesgue3D(N, r, s, t, 1000, 30);

return

function warp = Warpfactor1D(N, rout)

% Compute LGL and equidistant node distribution
LGLr = JacobiGL(0,0,N); req  = linspace(-1,1,N+1)';

% Compute V based on req
Veq = Vandermonde1D(N,req);

% Evaluate Lagrange polynomial at rout
Nr = length(rout); Pmat = zeros(N+1,Nr);
for i=1:N+1
    Pmat(i,:) = JacobiP(rout, 0, 0, i-1)';
end;
Lmat = Veq'\Pmat;

% Compute warp factor
warp = Lmat'*(LGLr - req);

return

function [warpx, warpy] = WarpShiftFace3D(p, pval1, pval2, pval3, L1,L2,L3,L4)

% function [warpx, warpy] = WarpShiftFace3D(p,pval, pval2, L1,L2,L3,L4)
% Purpose: compute warp factor used in creating 3D Warp & Blend nodes

[dtan1,dtan2] = evalshift(p, pval1, pval2, pval3, L2, L3, L4);
warpx = dtan1; warpy = dtan2;
return;


function [dx, dy] = evalshift(p, pval1, pval2, pval3, L1, L2, L3)

% function [dx, dy] = evalshift(p, pval, L1, L2, L3)
% Purpose: compute two-dimensional Warp & Blend transform

% 1) compute Gauss-Lobatto-Legendre node distribution
gaussX = -JacobiGL(0,0,p);

% 2) compute blending function at each node for each edge
blend1 = L2.*L3; blend2 = L1.*L3; blend3 = L1.*L2;

% 3) amount of warp for each node, for each edge
warpfactor1 = 4*evalwarp(p, gaussX, L3-L2);
warpfactor2 = 4*evalwarp(p, gaussX, L1-L3);
warpfactor3 = 4*evalwarp(p, gaussX, L2-L1);

% 4) combine blend & warp
warp1 = blend1.*warpfactor1.*(1 + (pval1*L1).^2);
warp2 = blend2.*warpfactor2.*(1 + (pval2*L2).^2);
warp3 = blend3.*warpfactor3.*(1 + (pval3*L3).^2);

% 5) evaluate shift in equilateral triangle
dx = 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3;
dy = 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;
return;
