function [u,time] = Advec2D(u, FinalTime)

Globals2D;
time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np,K); 

% compute time step size
rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

% outer time step loop 
while (time<FinalTime)
  
  if(time+dt>FinalTime), dt = FinalTime-time; end

   for INTRK = 1:5    
      % compute right hand side of TM-mode Maxwell's equations
      [rhsu] = AdvecRHS2D(u);

      % initiate and increment Runge-Kutta residuals
      resu = rk4a(INTRK)*resu + dt*rhsu;  
        
      % update fields
      u = u+rk4b(INTRK)*resu; 
   end;
   % Increment time
   time = time+dt;
end
return

function [rhsu] = AdvecRHS2D(u)

global vx
global vy

Globals2D;

% form field differences at faces
alpha=1; du = zeros(Nfp*Nfaces,K); 
du(:) = 0.5*(u(vmapM)-u(vmapP)).*(nx(:) - alpha*abs(nx(:)));
% 
% % impose boundary condition at x=0
% ubc  = exp(-1*( (Fx(mapB)-time).^2 + Fy(mapB).^2 + Fz(mapB).^2));
% du(mapB) = 0.5*(u(vmapB)-ubc).*(nx(mapB)-alpha*abs(nx(mapB)));
du(mapB) = 0;
% 
% compute right hand sides of the semi-discrete PDE
fx = vx.*u;
fy = vy.*u;
rhsu = -(rx.*(Dr*fx)+sx.*(Ds*fy)) + LIFT*(Fscale.*(du));

return;

