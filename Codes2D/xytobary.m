% convert from reference triangle coordinates to equilateral triangle
% coordinates 
%
% function [x, y] = rstoxy(r,s)
%
function [L1,L2,L3] = xytobary(x,y)

% equilateral vertices
v1 = 2*[-.5 -sqrt(3)/6]';
v2 = 2*[.5 -sqrt(3)/6]';
v3 = 2*[0 sqrt(3)/3]';

T = [v1(1)-v3(1) v2(1)-v3(1);
        v1(2)-v3(2) v2(2)-v3(2)];
for i = 1:length(x)    
    xy = [x(i);y(i)];
    L = T\[xy - v3];
    L1(i) = L(1);
    L2(i) = L(2);
    L3(i) = 1-L1(i)-L2(i);
end
L1 = L1(:);
L2 = L2(:);
L3 = L3(:);

return 