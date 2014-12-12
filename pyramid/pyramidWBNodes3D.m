% Interpolatory Warp and Blend nodes for the pyramid.
% 
% beta = optimizable blending parameter
function [r s t] = pyramidWBNodes3D(N,beta)
if nargin == 0
    N = 6;
    beta = 0;
end
if nargin ==1
    beta = 0; % optimization/blending parameter
end


% target nodes on surface - default to Nodes3D optimized blend for
% conformity with tet nodes
[mapr maps mapt] = pyramidSurfaceNodes3D(N);

% %% find coefficients mapping from equispaced to warped surface
[rbc sbc tbc] = pyramidSurfaceEquiNodes3D(N);
[Vbc v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, rbc, sbc, tbc, beta);

ids = [v_ids etri_ids equad_ids ftri_ids fquad_ids];
maprst = [mapr maps mapt];
mapcrst = Vbc(:,ids)\maprst;

% % evaluate map at equispaced volumed nodes
[req seq teq] = pyramidEquiNodes3D(N);
[Veq v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, req, seq, teq, beta);

ids = [v_ids etri_ids equad_ids ftri_ids fquad_ids];
rst = Veq(:,ids)*mapcrst; %% final coordinates of all nodes in triangle

r = rst(:,1); s = rst(:,2); t = rst(:,3);

%%
if nargin ==0
    plot3(mapr,maps,mapt,'ro');hold on
    plot3(r,s,t,'.');hold on;
    maxleb = PyramidLebesgue3D(N,r,s,t,5000,10);
    title(sprintf('leb const = %d',maxleb))
    view(-15,5)
end

function [V v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, r, s, t, beta)

%% pyramid
%%
%%    4------3
%%    |\    /|
%%    | \  / |
%%    |  5   |
%%    | /  \ |
%%    |/    \|
%%    1------2

% alphastore = [0;0;0;0.1002; 1.1332;1.5608;1.3413;1.2577;1.1603;...
%     1.10153;0.6080;0.4523;0.8856;0.8717;0.9655];
% alphaT = alphastore(N);

tol = 1e-10;
V = [];

% vertex functions
V(:,1) = .25*(1-r-s-t+r.*s./(1-t + tol));
V(:,2) = .25*(1+r-s-t-r.*s./(1-t + tol));
V(:,3) = .25*(1+r+s-t+r.*s./(1-t + tol));
V(:,4) = .25*(1-r+s-t-r.*s./(1-t + tol));
V(:,5) = t;

v_ids = 1:5;

sk = 6;

% edge function for triangular edges
edges = [1 5; 2 5; 3 5; 4 5];
etri_ids = [];
for e = 1:4
    i1 = edges(e,1);
    i2 = edges(e,2);
    
    switch e
        case 1
            blend1 = V(:,2);
            blend2 = V(:,4);
        case 2
            blend1 = V(:,1);
            blend2 = V(:,3);
        case 3
            blend1 = V(:,4);
            blend2 = V(:,2);
        case 4
            blend1 = V(:,1);
            blend2 = V(:,3);
    end
    blend2 = (1 + beta*(blend1).^2).*(1 + beta*(blend2).^2); % blend towards opposite vertex on faces
    
    for i=0:N-2
        xi = V(:,i1)-V(:,i2);
        V(:,sk) = blend2.*V(:,i1).*V(:,i2).*JacobiP(xi, 1, 1, i);
        etri_ids = [etri_ids sk];
        sk = sk+1;
    end
end

% edge functions for base
edges = [1 2; 2 3; 3 4; 4 1];
equad_ids = [];
for e = 1:4
    i1 = edges(e,1);
    i2 = edges(e,2);
    
    %     op_vert = setdiff(1:5,edges(e,:));
    %V(:,op_vert(1)).*V(:,op_vert(2)).*V(:,op_vert(3));
    blend2 = (1 + beta*(V(:,5)).^2); % blend towards top vertex
    
    for i=0:N-2
        xi = V(:,i1)-V(:,i2);
        V(:,sk) = blend2.*V(:,i1).*V(:,i2).*JacobiP(xi, 1, 1, i);
        equad_ids = [equad_ids sk];
        sk = sk+1;
    end
end

%triangular faces
faces = [1 2 5; 2 3 5; 3 4 5; 4 1 5];
ftri_ids = [];
for f = 1:4
    i1 = faces(f,1);
    i2 = faces(f,2);
    i3 = faces(f,3);
    
    % bubble edge blend
    %         op_vert = setdiff(1:5,faces(f,:));  blend2 = V(:,op_vert(1)).*V(:,op_vert(2));
    
    % plane edge blend
    switch f;
        case 1;
            blend2 = ((1+s)-t)/2;
        case 2;
            blend2 = ((1-r)-t)/2;
        case 3;
            blend2 = ((1-s)-t)/2;
        case 4;
            blend2 = ((1+r)-t)/2;
    end
    %blend2 = 1-V(:,5); % blend to base
    
    blend2 = 1 + beta*(blend2).^2;
    %blend2 = beta*(V(:,5)-.5); % blend to base - switches negative at midway pt.
    
    
    L1 = V(:,i1);    L2 = V(:,i2);    L3 = V(:,i3);
    [x y] = eqbarytoxy(L1,L2,L3);
    [rr ss] = xytors(x,y);
    Vf = Vandermonde2D(N-3,rr,ss);
    for i = 1:size(Vf,2)
        V(:,sk) = blend2.*V(:,i1).*V(:,i2).*V(:,i3).*Vf(:,i);
        ftri_ids = [ftri_ids sk];
        sk = sk + 1;
    end
end

% square face on the bottom
blend2 = 1 + beta*(V(:,5)).^2;
fquad_ids = [];
for i = 0:N-2
    for j = 0:N-2
        V(:,sk) = blend2.*V(:,1).*V(:,2).*V(:,3).*V(:,4).*JacobiP(r,1,1,i).*JacobiP(s,1,1,j);
        fquad_ids = [fquad_ids sk];
        sk = sk + 1;
    end
end


return;



