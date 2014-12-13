% convert from reference triangle coordinates to equilateral triangle
% coordinates 
%
% function [x, y] = rstoxy(r,s)
%
function [L1,L2,L3,L4] = rsttobary(rin,sin,tin)

% rst vertices
[x y z] = EquiNodes3D(1);
v1 = [x(1);y(1);z(1)];
v2 = [x(2);y(2);z(2)];
v3 = [x(3);y(3);z(3)];
v4 = [x(4);y(4);z(4)];

T = [x(1)-x(4), x(2)-x(4), x(3)-x(4);
     y(1)-y(4), y(2)-y(4), y(3)-y(4);
     z(1)-z(4), z(2)-z(4), z(3)-z(4)];
 
for i = 1:length(rin)
    xy = [rin(i);sin(i);tin(i)];
    L = T\[xy - v4];
    L1(i) = L(1);
    L2(i) = L(2);
    L3(i) = L(3);
    L4(i) = 1-L1(i)-L2(i)-L3(i);
end
L1 = L1(:);
L2 = L2(:);
L3 = L3(:);
L4 = L4(:);

return 