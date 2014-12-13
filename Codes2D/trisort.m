% sort triangle xy coordinates from left to right, bottom to top

function [xsorted ysorted] = trisort(x,y,N)

xsorted=[];
ysorted=[];
n = N+1;
for i = 1:N+1 % number of rows 
    [d py] = sort(y);
    py = py(1:n);
    [d px] = sort(x(py));    
    py = py(px);    
    xsorted = [xsorted; x(py)];
    ysorted = [ysorted; y(py)];
    newInds = setdiff(1:length(x),py);
    x = x(newInds);
    y = y(newInds);
    n=n-1; % number of nodes in a row
end