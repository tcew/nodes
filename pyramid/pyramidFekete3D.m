% ODE-based steepest ascent for the Fekete nodes.  
% 
% http://www.c3.lanl.gov/~wingate/papers.html
% http://www.scd.ucar.edu/css/staff/taylor/points/fekete.html

function [r,s,t] = pyramidFekete3D(inN, rin, sin, tin, tol)

Globals3D

if nargin ==0
    inN = 10;
    %     [r s t] = pyramidEquiNodes(inN);
    [rin sin tin] = pyramidWBNodes3D(inN); % initialize to Warp and Blend nodes
end
if nargin < 4
    [rin sin tin] = pyramidWBNodes3D(inN); % initialize to Warp and Blend nodes
end
if nargin < 5
    tol = 1e-10;
end

% set the nodes to be equispaced
N = inN; r = rin; s = sin; t = tin;

Np = numel(r);

% location of vertices and edges
% tol = 1e-10;
bmask = abs(t) < tol; % bottom boundary
bmask = bmask + (abs(t + r - 1) < tol | abs(t - r - 1) < tol);
bmask = bmask + (abs(t + s - 1) < tol | abs(t - s - 1) < tol);
bmask = bmask > 0; % boolean

% storage for node velocities
U = zeros(Np,1);
V = zeros(Np,1);
W = zeros(Np,1);

% RK storage for X and Zeta time integrating residuals
resR = zeros(Np,1);
resS  = zeros(Np,1);
resT  = zeros(Np,1);

% constants
%RKtime = 0;
time   = 0;
dt     = 0.5/(N^3);
stoptol = tol/dt
FinalTime = 1000;
Nsteps = FinalTime/dt;

% outer time step loop
for tstep = 1:Nsteps
    
    % inner multi-stage Runge-Kutta loop
    for INTRK = 1:5
        % initiate Runge-Kutta stage
        resR   = rk4a(INTRK)*resR;
        resS   = rk4a(INTRK)*resS;
        resT   = rk4a(INTRK)*resT;
        %RKtime = time+dt*rk4c(INTRK);
        
        % calculate derivative matrices
        [VDM Dr Ds Dt] = pyramidBasisBergot3D(N, r, s, t);
        if tstep == 1 && INTRK ==1
            if nargin==0
                maxleb = PyramidLebesgue3D(N, r, s, t, 20000, 15);
                disp(sprintf('det at t=0: %d, leb const = %d',det(VDM),maxleb))
            else
                disp(sprintf('det at t=0: %d',det(VDM)))
            end            
        end
        invV = inv(VDM);
        
        Dr = Dr * invV;
        Ds = Ds * invV;
        Dt = Dt * invV;
        
        % set velocity of each node, to the gradient of
        % its Lagrange interpolating polynomial
        U = diag(Dr);
        V = diag(Ds);
        W = diag(Dt);
        
        % constrain boundary nodes
        U(bmask) = 0;
        V(bmask) = 0;
        W(bmask) = 0;
        
        % update R-K residuals
        resR = resR+dt*U;
        resS = resS+dt*V;
        resT = resT+dt*W;
        
        % finish Runge-Kutta stage
        r = r + rk4b(INTRK)*resR;
        s = s + rk4b(INTRK)*resS;
        t = t + rk4b(INTRK)*resT;
    end
    % end inner-RK loop
    
    if mod(tstep,100)==0
        disp(sprintf('tstep %i out of Nsteps = %i, maxU = %d, maxV = %d, maxW = %d',tstep,Nsteps,maxU,maxV,maxW))            
    end
    
    if nargout == 0
        if(tstep == 0 || mod(tstep,10) == 0)
            clf
            color_line3(r(bmask),s(bmask),t(bmask),t(bmask),'o');
            color_line3(r(~bmask),s(~bmask),t(~bmask),t(~bmask),'.');
            hold off;
            title(sprintf('maxU = %d, maxV = %d, maxW = %d',maxU,maxV,maxW))
            view(3)
            axis([-1 1 -1 1 0 1])
            
            drawnow;
        end
    end
    % update time
    time = time+dt;
    
    % stopping criterion
    maxU = max(abs(U));
    maxV = max(abs(V));
    maxW = max(abs(W));
    
    if ((maxU < stoptol) && (maxV < stoptol) && (maxW < stoptol))
        break;
    end
end
VDM = pyramidBasisBergot3D(N, r, s, t);
if nargin ==0
    maxleb = PyramidLebesgue3D(N, r, s, t, 20000, 10);
    %disp(sprintf('det at tfinal: %d, leb const = %d',det(VDM),maxleb))
    keyboard
else
    disp(sprintf('det at tfinal: %d',det(VDM)))
end




