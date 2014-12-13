% convert from reference triangle coordinates to equilateral triangle
% coordinates 
%
% function [x, y] = rstoxy(r,s)
%
function [L1,L2,L3,L4] = xyztobary(xin,yin,zin)

% equilateral vertices
v1 = [-1, -1/sqrt(3), -1/sqrt(6)]; v2 = [ 1, -1/sqrt(3),-1/sqrt(6)];
v3 = [ 0,  2/sqrt(3), -1/sqrt(6)]; v4 = [ 0;  0;         3/sqrt(6)];

x = [v1(1) v2(1) v3(1) v4(1)];
y = [v1(2) v2(2) v3(2) v4(2)];
z = [v1(3) v2(3) v3(3) v4(3)];

T = [x(1)-x(4), x(2)-x(4), x(3)-x(4);
     y(1)-y(4), y(2)-y(4), y(3)-y(4);
     z(1)-z(4), z(2)-z(4), z(3)-z(4)];
 
for i = 1:length(xin)
    xy = [xin(i);yin(i);zin(i)];
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