function D = spdiag(d)
n = length(d(:));
D = spdiags(d(:),0,n,n);

