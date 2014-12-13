% convert from reference triangle coordinates to equilateral triangle
% coordinates 
%
% function [x, y] = rstoxy(r,s)
%
function [x,y] = eqbarytoxy(L1,L2,L3)

% equilateral vertices
v1 = 2*[-.5 -sqrt(3)/6]';
v2 = 2*[.5 -sqrt(3)/6]';
v3 = 2*[0 sqrt(3)/3]';

for i = 1:length(L1)
    XY = v1*L1(i) + v2*L2(i) + v3*L3(i);
    x(i) = XY(1);
    y(i) = XY(2);
end

return 