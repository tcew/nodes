% uses Sommariva/Vianello's greedy optimization algorithm for discrete Fekete
% points on the pyramid

function [r s t] = pyramidApproxFekete_greedy(N,Nsample)

% N = 5;
if nargin==1
    Nsample = N^2+1;
end
[r s t] = pyramidEquiNodes3D(Nsample);

% [rq sq tq w] = pyramidCubature3D(N+2); % for normalization
%Vq = pyramidBasisBergot3D(N,rq,sq,tq);
%invNorms = spdiag(1./sqrt(diag(Vq'*spdiag(w)*Vq)));

tol = 1e-10;

[rb sb tb] = pyramidSurfaceNodes3D(N);
Nb = length(rb);

nlevels = 1;%15;
radius = norm([r(1) s(1) t(1)] - [r(2) s(2) t(2)]);

for level = 1:nlevels
    % skim off boundary nodes
    bmask = abs(t) < tol; % bottom boundary
    bmask = bmask + (abs(t + r - 1) < tol | abs(t - r - 1) < tol);
    bmask = bmask + (abs(t + s - 1) < tol | abs(t - s - 1) < tol);
    bmask = bmask > 0; % boolean
    r(bmask) = []; s(bmask) = []; t(bmask) = [];
    r = [rb(:);r(:)];
    s = [sb(:);s(:)];
    t = [tb(:);t(:)];
           
    % VDM matrix
    V = pyramidBasisBergot3D(N,r,s,t); 
    for j = 1:size(V,2)
        V(:,j) = V(:,j)/norm(V(:,j)); %l2 normalize cols
    end
    
    % get ids according to greedy algo
    ids = Greedy(V',Nb);
    
    if nargout == 0
        maxleb = PyramidLebesgue3D(N, r(ids), s(ids), t(ids), 10000, 15);
        disp(sprintf('leb constant at level %i is %d',level,maxleb))
        
        plot3(r(ids),s(ids),t(ids),'.');
        title(sprintf('level = %i, leb const = %d',level,maxleb));
        drawnow
    end
    
    if level < nlevels
        Nsamp = 1000;
        %rad = linspace(-radius,radius,Nsamp);
        %[theta, phi] = meshgrid(linspace(0,2*pi,Nsamp));
%         rad = radius*randn(Nsamp);
%         [theta, phi] = meshgrid(2*pi*rand(Nsamp),2*pi*rand(Nsamp));
%         
%         rn = rad(:)*(cos(theta(:)).*sin(phi(:)))';
%         sn = rad(:)*(sin(theta(:)).*sin(phi(:)))';
%         tn = rad(:)*cos(phi(:))';
%         rn = rn(:); sn = sn(:); tn = tn(:);
        rn = radius*randn(Nsamp,1);
        sn = radius*randn(Nsamp,1);
        tn = radius*randn(Nsamp,1);
        
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
        
        radius = radius/2        
%         keyboard
    end   
    size(r)
end
r = r(ids);
s = s(ids);
t = t(ids);

function piv = Greedy(A,K)

% assume N<=M
N = size(A,1);
M = size(A,2);

% use first K columns
piv = 1:K;
for k=1:K
    col = A(:,k)/norm(A(:,k));
    
    % orthogonalize
    A = A - col*(col'*A);    
end

for k=K+1:N
    normA = sum(A.^2,1);
    [foo,ik] = max(normA);
    
    col = A(:,ik)/norm(A(:,ik));
    
    % orthogonalize
    A = A - col*(col'*A);
    
    piv = [piv,ik];
    k/N;
end
piv = unique(piv);
