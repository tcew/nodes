% nodeSwitch = 1: W&B
% nodeSwitch = 2: GQ
% nodeSwitch = 3: Williams sym pts
%
function AdvectionDriver2D(nodeSwitch,Nin,FinalTime,filename)

global plotFlag

% plot if save data
plotFlag = 0;
saveData=1;
if (nargin<4)
    saveData=0;
    plotFlag=1;
end
% Driver script for solving the 2D advection
Globals2D;

% Polynomial order used for approximation 
% N = 6
N = Nin % must do this b/c of issues w/global variables 

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');

% Initialize solver and construct grid and metric
StartUp2D;

%rebuild maps
BuildPeriodicMaps2D(2,2);

% plot(x,y,'.')
% hold on
% for i = 1:length(vmapP)
%     plot(x(vmapM(i)),y(vmapM(i)),'ro')
%     pause
%     plot(x(vmapP(i)),y(vmapP(i)),'ro')
%     pause
% end
% keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Interp
global InterpBack
global xq 
global yq
global F 

% compute rq, sq based on nodal points
if (nodeSwitch==1) % warp and blend 
    rq = r;
    sq = s;    
elseif (nodeSwitch==2) % get GQ points (if exist)
    NcubOrder = findNcubOrder(N);
    [rq,sq,wq, Ncub] = Cubature2D(NcubOrder);
    error('No quad rule matching')    
elseif (nodeSwitch==3) % Williams symmetric points    
    [xwq,ywq,wwq] = QNodes2D(N); 
    [rq,sq] = xytors(xwq,ywq);
elseif (nodeSwitch==4) % experimental nodes    
    [xjq,yjq] = blyth(N,JacobiGQ(0,0,N));
    [rq,sq] = xytors(xjq,yjq);
end

% make interp mats and qpoints
xq = 0.5*(-(rq+sq)*VX(va)+(1+rq)*VX(vb)+(1+sq)*VX(vc));
yq = 0.5*(-(rq+sq)*VY(va)+(1+rq)*VY(vb)+(1+sq)*VY(vc));
Vq = Vandermonde2D(N,rq,sq); invVq = inv(Vq);
Interp = Vq*invV; % interp to "q" nodes
InterpBack = V*invVq; % interp from "q" nodes to WaB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 EXTRA QUADRATURE  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global useCubature;
useCubature = 1
if (useCubature)
    cubDeg = findNcubOrder(N);
    [cubR,cubS,cubW, Ncub] = Cubature2D(cubDeg);
%     [cubX,cubY,cubW] = QNodes2D(N+1);[cubR cubS] = xytors(cubX,cubY);
%     [cubX,cubY,cubW] = QNodes2D(N,cubDeg);[cubR cubS] = xytors(cubX,cubY);
    Ncub = length(cubW);
    
    disp(['Np = ', num2str(Np), ', Ncub = ', num2str(Ncub)])
       
    Vcub = Vandermonde2D(N,cubR,cubS);    
    Interp = Vcub*invV; % interp to cubature points
    Wcub = diag(cubW);    
    Mcub = Interp'*Wcub*Interp; invMcub = inv(Mcub);
    InterpBack = invMcub*Interp'*Wcub; % InterpBack*diag(vx(:))*Interp does a projection
    er = norm(eye(Np)-InterpBack*Interp,'fro');    
    disp(['integration err in mass mat = ' num2str(er)])
    
    xq = 0.5*(-(cubR+cubS)*VX(va)+(1+cubR)*VX(vb)+(1+cubS)*VX(vc));
    yq = 0.5*(-(cubR+cubS)*VY(va)+(1+cubR)*VY(vb)+(1+cubS)*VY(vc));    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global vx
global vy

% vx = @(x,y) ones(size(x));

% vx = @(x,y) .1 + (1-y.^2).^5;
vx = @(x,y) .1 + exp(-10*y.^2);
vy = @(x,y) zeros(size(x));

% a = .5;
% vx = @(x,y) (1-a)*(.1 + (1-y.^2).^5);
% vy = @(x,y) a*(.1 + (1-x.^2).^5);

% vx = @(x,y) -cos(pi/2*x).*sin(pi/2*y); 
% vy = @(x,y) sin(pi/2*x).*cos(pi/2*y); 
% vx = @(x,y) -y;
% vy = @(x,y) x;

t = pi/4;
o = @(x,y) -sin(t)*x + cos(t)*y;
b = @(x,y) (1-x.^2).^2.*(1-y.^2).^2; % bubble 
c = .1;
% vx = @(x,y) cos(t)*(c + b(x,y).*exp(-10*o(x,y).^2))/(exp(0)+c);
% vy = @(x,y) sin(t)*(c + b(x,y).*exp(-10*o(x,y).^2))/(exp(0)+c);

vx = vx(xq,yq);
vy = vy(xq,yq);

% filtering
F = eye(Np);
% F = Filter2D(N,N-2,.1); % b/w .01-.25
% F = CutOffFilter2D(N-1,N/(2*N+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Set initial conditions
cx = .1;
cy = .1;
D = (x-cx).^2 + (y-cy).^2;
u = exp(-D*5^2).*(1-x.^2).*(1-y.^2);
% u = ones(size(x));

global alpha
alpha = 1;

if (Np*K<5000)
    u0 = u;
    [RHSDG rhsVol rhsSurf] = getOperator(u0);
    rhsuOp = RHSDG*u0(:);
    [rhsu0] = AdvecRHS2D(u0,0);
    disp(['norm of diff = ', num2str(norm(rhsuOp-rhsu0(:)))])
end
if (size(RHSDG,1)<5000)
    max(real(eigs(RHSDG,100,'lr')))
elseif (size(RHSDG,1)<1000)
    max(real(eig(full(RHSDG))))
else
    disp(['matrix size = ',num2str(size(RHSDG,1))])
end

keyboard

% return
% u = (1-y.^2).*sin(pi*x);
% PlotField2D(N, x, y, u); 
% colorbar
% view(2)

% 
% % Solve Problem
% FinalTime = 125;

% figure(1); 
[u,time,uNorm] = Advec2D(u,FinalTime);
% keyboard
% figure
% plot(uNorm)
% PlotField2D(N, x, y, u); 
% colorbar
% view(2)
if (saveData)    
    save(filename,'time','uNorm')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NcubOrder = findNcubOrder(N)
Np = (N+1)*(N+2)/2;
NcubOrder = 0;
for i = 1:28
    [cubR,cubS,cubW, Ncub] = Cubature2D(i);
    if (Ncub>Np) % find set w/more points than nodal pts
        NcubOrder = i;
        return;
    end
end

function [RHSDG rhsVol rhsSurf] = getOperator(u)

Globals2D

global vx
global vy
global Interp;
global InterpBack;

blkInterp = kron(speye(K),Interp);
blkInterpBack = kron(speye(K),InterpBack);
NqC = size(blkInterp,1);

blkFx = blkInterpBack*spdiags(vx(:),0,NqC,NqC)*blkInterp;
blkFy = blkInterpBack*spdiags(vy(:),0,NqC,NqC)*blkInterp;

blkDr = kron(speye(K),Dr);
blkDs = kron(speye(K),Ds);
global F
blkFilter = kron(speye(K),F);
blkLIFT = kron(speye(K),LIFT);

% extraction ops for boundary dofs - Em*u(:) = u(vmapM)
NfpK = Nfp*Nfaces*K;
Em = sparse(1:NfpK,vmapM,ones(NfpK,1),NfpK,Np*K); 
Ep = sparse(1:NfpK,vmapP,ones(NfpK,1),NfpK,Np*K); 

% divu = rx.*Dr*Fx + sx.*Ds*Fx + ry.*Dr*Fy + sy.*Ds*Fy;
% dudx = rx.*Dr*u + sx.*Ds*u , dudy = ry.*Dr*u + sy.*Ds*u 
Dx = spdiags(rx(:),0,Np*K,Np*K)*blkDr + spdiags(sx(:),0,Np*K,Np*K)*blkDs;
Dy = spdiags(ry(:),0,Np*K,Np*K)*blkDr + spdiags(sy(:),0,Np*K,Np*K)*blkDs;
rhsVol = -Dx*blkFx - Dy*blkFy;
% rhsVol = -blkFx*Dx - blkFy*Dy; % beta*grad(u) for div-free beta

global alpha
alphax = (nx(:)-(1-alpha)*abs(nx(:)))/2;
alphay = (ny(:)-(1-alpha)*abs(ny(:)))/2;

dU = spdiags(alphax,0,NfpK,NfpK)*(Em-Ep)*blkFx + ...
    spdiags(alphay,0,NfpK,NfpK)*(Em-Ep)*blkFy;
rhsSurf = blkLIFT*spdiags(Fscale(:),0,NfpK,NfpK)*dU;

RHSDG = blkFilter*(rhsVol+rhsSurf);
% RHSDG(vmapB,vmapB) = eye(length(vmapB(:)));

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,time,uNorm] = Advec2D(u, FinalTime)

global plotFlag

Globals2D;
time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np,K); 

% compute time step size
rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

% outer time step loop 
i = 0;
while (time<FinalTime)
  
  if(time+dt>FinalTime), dt = FinalTime-time; end

   for INTRK = 1:5    
      % compute right hand side of TM-mode Maxwell's equations
      [rhsu] = AdvecRHS2D(u,time);

      % initiate and increment Runge-Kutta residuals
      resu = rk4a(INTRK)*resu + dt*rhsu;  
        
      % update fields
      u = u+rk4b(INTRK)*resu; 
   end;
   % Increment time
   time = time+dt;
   i = i+1;   
   Mu = MassMatrix*(J.*u);
   uNorm(i) = sqrt(u(:)'*Mu(:));   
   if (mod(i,25)==0)
       disp(['t = ', num2str(time), ', norm of u = ', num2str(uNorm(i))])
       if plotFlag
           PlotField2D(N, x, y, u);
           colorbar
           caxis([-.1,1.1])
           view(2)
           title(['Time = ', num2str(time)])
           drawnow
       end
   end
end

return

function [rhsu] = AdvecRHS2D(u,t)

global vx
global vy

global Interp
global InterpBack
global xq
global yq

global F

Globals2D;

fx = InterpBack*(vx.*(Interp*u));
fy = InterpBack*(vy.*(Interp*u));

% form field differences at faces
% alpha = 1 = central, 0 = upwind
fxF = fx(vmapM)-fx(vmapP);
fyF = fy(vmapM)-fy(vmapP);

global alpha
du = zeros(Nfp*Nfaces,K); 
du(:) = .5*((nx(:)-(1-alpha)*abs(nx(:))).*fxF + (ny(:)-(1-alpha)*abs(ny(:))).*fyF);
% du(mapB) = 0;

% 
% compute right hand sides of the semi-discrete PDE
rhsu = -Div2D(fx,fy) + LIFT*(Fscale.*(du));
rhsu = F*rhsu;
% beta*grad(u) version for div-free discretizations
% [dudx dudy] = Grad2D(u); 
% rhsu = -InterpBack*(vx.*(Interp*dudx) + vy.*(Interp*dudy)) + LIFT*(Fscale.*(du));

return
