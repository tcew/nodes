% % Driver script for solving the 2D advection
% Globals2D;
% % Polynomial order used for approximation
% N = 6
% % Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% % Initialize solver and construct grid and metric
% StartUp2D;

x = linspace(-1,1,2048);
[x y] = meshgrid(x);

vx = @(x,y,t) (.1 + (1-y.^2).^5);

vx = @(x,y,t) (1 + .5*sin(pi*t)*(1-y.^2).^5);

vy = @(x,y) zeros(size(x));

figure
for t = 0:.1:5
    xt = x - vx(x,y,t).*t; % stretched xcoord
    xtp = mod(xt+1,2)-1; % periodicity
    D = (xtp).^2 + y.^2;
    u = exp(-D*5^2);
    %     PlotField2D(N, x, y, u);
%     surf(x,y,u)
plot(xt,y)
    axis equal;shading interp;view(2)
    title(['time = ', num2str(t)])
    drawnow    
end
