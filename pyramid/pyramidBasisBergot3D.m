% orthonormalized rational pyramid basis of Bergot, Durufle, Morgane.

function [V, Dr, Ds, Dt] = pyramidBasisBergot3D(N,r,s,t)

[V Dr Ds Dt] = bergotBasis(N,r,s,t);

function [V Dr Ds Dt] = bergotBasis(N,r,s,t)

tol = 1e-12;
t(abs(t-1)<tol) = 1-tol;

Nfp = 3*N^2+2;
Nip = 1/6*(N-1)*(N-2)*(2*N-3);
Np = Nfp + Nip;

% Vandermonde matrix
V = zeros(numel(r),Np);
if nargout > 1
    Dr = V; Ds = V; Dt = V;
end

% precompute for speed 
p1 = cell(N,1); dp1dr = cell(N,1);
p2 = cell(N,1); dp2ds = cell(N,1);
jp3 = cell(N,N); gjp3 = cell(N,N);
for i = 0:N
    p1{i+1} = JacobiP(r./(1-t),0,0,i);
    dp1dr{i+1} = GradJacobiP(r./(1-t + tol), 0, 0, i);
    
    % at t = 1, r,s = 0
    p1{i+1}(abs(t-1)<tol) = JacobiP(0,0,0,i);
    dp1dr{i+1}(abs(t-1)<tol) = GradJacobiP(0, 0, 0, i);
    
    p2{i+1} = JacobiP(s./(1-t),0,0,i);
    dp2ds{i+1} = GradJacobiP(s./(1-t), 0, 0, i);
    
    % at t = 1, r,s = 0
    p2{i+1}(abs(t-1)<tol) = JacobiP(0,0,0,i);
    dp2ds{i+1}(abs(t-1)<tol) = GradJacobiP(0, 0, 0, i);

    for mij = 0:N
        jp3{i+1}{mij+1} = JacobiP(2*t-1,2*mij + 2,0,i);
        gjp3{i+1}{mij+1} = GradJacobiP(2*t-1, 2*mij+2, 0, i);
    end
end

off = 1;
for i = 0:N    
    for j = 0:N        
        mij = max(i,j);
        for k = 0:N-mij
            p3 = (1-t).^mij.*jp3{k+1}{mij+1};
            dp3dt = -mij*(1-t).^(mij-1).*jp3{k+1}{mij+1} + (1-t).^mij.*2.*gjp3{k+1}{mij+1};                        
            
            % in case (1-t) = 0.  
            if mij ==0
                dp3dt = (1-t).^mij.*2.*gjp3{k+1}{mij+1};
            end
            
            scale = sqrt(2*4.^(mij+1));
            V(:,off) = p1{i+1}.*p2{j+1}.*p3*scale;
            if nargout > 1
                Dr(:,off) = dp1dr{i+1}.*p2{j+1}.*p3;
                Ds(:,off) = p1{i+1}.*dp2ds{j+1}.*p3;
                Dt(:,off) = p1{i+1}.*p2{j+1}.*dp3dt;
            end
            
            off = off + 1;
        end
    end
end
