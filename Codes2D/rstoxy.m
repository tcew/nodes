% convert from reference triangle coordinates to equilateral triangle
% coordinates 
%
% function [x, y] = rstoxy(r,s)
%
function [x, y] = rstoxy(r,s)

% reference triangle coordinates
v1 = [-1 -1]';
v2 = [1 -1]';
v3 = [-1 1]';

% T = [v1(1)-v3(1) v2(1)-v3(1);
%      v1(2)-v3(2) v2(2)-v3(2)];
% for i = 1:length(r)    
%     rs = [r(i);s(i)];
%     L = T\[rs - v3];
%     L1(i) = L(1);
%     L2(i) = L(2);
%     L3(i) = 1-L1(i)-L2(i);
% end

denom = (v2(2)-v3(2))*(v1(1)-v3(1)) + (v3(1)-v2(1))*(v1(2)-v3(2));
L1 = ((v2(2)-v3(2))*(r-v3(1)) + (v3(1)-v2(1))*(s-v3(2)))/denom;
L2 = ((v3(2)-v1(2))*(r-v3(1)) + (v1(1)-v3(1))*(s-v3(2)))/denom;
L3 = 1-L1-L2;

% equilateral vertices
v1 = 2*[-.5 -sqrt(3)/6]';
v2 = 2*[.5 -sqrt(3)/6]';
v3 = 2*[0 sqrt(3)/3]';

for i = 1:length(r)
    XY = v1*L1(i) + v2*L2(i) + v3*L3(i);
    x(i,1) = XY(1);
    y(i,1) = XY(2);
end

return 