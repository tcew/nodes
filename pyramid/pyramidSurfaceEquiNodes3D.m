% constructs equispaced nodes on the pyramid surface. 

function [x,y,z] = pyramidSurfaceEquiNodes3D(N)


% build surface nodes for a degree 5 pyramid trace
% [-1,1]x[-1,1]x[0,1] with top vertex at (0,0,1)

% start with W&B nodes for triangle
%[r,s] = Nodes2D(N);
%[r,s] = xytors(r,s);
[r s t] = EquiNodes3D(N);

inds = abs(t+1)<1e-8;
r = r(inds);s = s(inds);

% r = r(:);s = s(:);
Np = length(r);

tol = 1e-8;
enums(:,1) = find(abs(s+1)<tol);
enums(:,2) = find(abs(s+r)<tol);
enums(:,3) = find(abs(r+1)<tol);

eids = unique(enums(:));
iids = setdiff(1:Np, eids);

% 
[gr] = linspace(-1,1,N+1);

[r2d,s2d] = meshgrid(gr);
t2d = zeros(N+1,N+1);

r2d = r2d(:);
s2d = s2d(:);
t2d = t2d(:);

I = ones(Np,1);

%% Face 2
x1 = -1; y1 = -1; z1 = +0;
x2 =  1; y2 = -1; z2 = +0;
x3 =  0; y3 =  0; z3 = +1;
x3d = [r2d;0.5*(-(r+s)*x1 + (1+r)*x2 + (1+s)*x3)];
y3d = [s2d;0.5*(-(r+s)*y1 + (1+r)*y2 + (1+s)*y3)];
z3d = [t2d;0.5*(-(r+s)*z1 + (1+r)*z2 + (1+s)*z3)];


%% Face 3
x1 =  1; y1 = -1; z1 = +0;
x2 =  1; y2 = +1; z2 = +0;
x3 =  0; y3 =  0; z3 = +1;
x3d = [x3d;0.5*(-(r+s)*x1 + (1+r)*x2 + (1+s)*x3)];
y3d = [y3d;0.5*(-(r+s)*y1 + (1+r)*y2 + (1+s)*y3)];
z3d = [z3d;0.5*(-(r+s)*z1 + (1+r)*z2 + (1+s)*z3)];

%% Face 4
x1 =  1; y1 = +1; z1 = +0;
x2 = -1; y2 = +1; z2 = +0;
x3 =  0; y3 =  0; z3 = +1;
x3d = [x3d;0.5*(-(r+s)*x1 + (1+r)*x2 + (1+s)*x3)];
y3d = [y3d;0.5*(-(r+s)*y1 + (1+r)*y2 + (1+s)*y3)];
z3d = [z3d;0.5*(-(r+s)*z1 + (1+r)*z2 + (1+s)*z3)];
 
%% Face 5
x1 =  -1; y1 = -1; z1 = +0;
x2 =  -1; y2 = +1; z2 = +0;
x3 =  0; y3 =  0; z3 = +1;
x3d = [x3d;0.5*(-(r+s)*x1 + (1+r)*x2 + (1+s)*x3)];
y3d = [y3d;0.5*(-(r+s)*y1 + (1+r)*y2 + (1+s)*y3)];
z3d = [z3d;0.5*(-(r+s)*z1 + (1+r)*z2 + (1+s)*z3)];

%% remove duplicates
tol = 1e-10;
cnt = 1;
x = x3d(1); y = y3d(1); z = z3d(1);
for n=2:length(x3d)
    d = (x3d(n)-x).^2 + (y3d(n)-y).^2 + (z3d(n)-z).^2;
    
    if(min(d)>tol)
        x = [x;x3d(n)];
        y = [y;y3d(n)];
        z = [z;z3d(n)];
    end
end

% scatter3(x,y,z)

pyrNp = 5 + (N-1)*8 + 4*(N-2)*(N-1)/2 + (N-1)*(N-1);
calcPyrNp = length(x);



