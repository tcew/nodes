% uses Sommariva/Vianello's iterative refinement algorithm for discrete Fekete
% points on the pyramid

function [r s t] = pyramidApproxFekete_QR(N,Nsample)

tol = 1e-8;

[rb sb tb] = pyramidSurfaceNodes3D(N);

if N < 3
    r = rb; s = sb; t = tb;
    return
end

if nargin == 1
    Nsample = N^2+1;
end

[r s t] = pyramidEquiNodes3D(Nsample);

[rq sq tq wq] = pyramidCubature3D(N+2);
Vq = pyramidBasisBergot3D(N,rq,sq,tq);

% remove boundary part of VDM
VDMb = pyramidBasisBergot3D(N,rb,sb,tb);
[Ub,Sb,Vb] = svd(VDMb);
Nb = numel(rb);
Vbb = Vb(:,1:Nb); Sbb = Sb(:,1:Nb); Vbi = Vb(:,Nb+1:end);

nlevels = 1;
radius = norm([r(1) s(1) t(1)] - [r(2) s(2) t(2)]);
for level = 1:nlevels              
    
    % skim off boundary nodes
    bmask = abs(t) < tol; % bottom boundary
    bmask = bmask + (abs(t + r - 1) < tol | abs(t - r - 1) < tol);
    bmask = bmask + (abs(t + s - 1) < tol | abs(t - s - 1) < tol);
    bmask = bmask > 0; % boolean
    r(bmask) = []; s(bmask) = []; t(bmask) = [];
    
    % take interior modes
    V = pyramidBasisBergot3D(N,r,s,t);
    V = V*Vbi; 
    
    P = eye(size(V,2));    
    
    niter = 4;
    for k = 0:niter
        [~, R] = qr(V,0);
        U = inv(R);
        V = V*U;
        P = P*U;
    end
        
    m = sum((Vq*Vbi)'*spdiag(wq),2);
        
    mu = (P')\m;
    w = (V')\mu;
    ids = find(abs(w)>tol);
    
    if nargin ==0
        rp = [rb(:); r(ids(:))];
        sp = [sb(:); s(ids(:))];
        tp = [tb(:); t(ids(:))];
        plot3(rp,sp,tp,'.');drawnow
    end
    
    if level < nlevels
        Nsamp = 10;
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
        
        radius = radius/4        
    end        
    
end

r = [rb(:); r(ids(:))];
s = [sb(:); s(ids(:))];
t = [tb(:); t(ids(:))];
