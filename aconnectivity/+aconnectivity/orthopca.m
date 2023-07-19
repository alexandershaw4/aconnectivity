function [U,V] = orthopca(x,k,opt)
% factorize matrix x(m*n) into k components, returning
% U(mxk) and V(kxn) where each V(ki,n) can only be non-zero for one component (k)
%
% usage: [U,V] = orthopca(x,k)
%
% AS2023

if nargin < 3 || isempty(opt)
    opt = 'svd';
end

% select initial factorisation algorithm
switch opt ...
    case 'svd';

        [u,s,v] = svd(x,'econ');
        
        U = u(:,1:k);
        S = s(1:k,1:k);
        V = v(:,1:k);
        U = U*S;
    case 'nnmf'
        [U,V] = nnmf(x,k);
        V = V';
end

% max number iterations
iter = 64;

for i = 1:iter;
    
    % start
    UX = U; VX = V;

    % single component membership
    V  = singlemember(V);
    r0 = norm(x - (U*V'));

    % update u
    b = (V\x')';

    % update v
    bv = x'/b';
    bv = bv.*~~V;
    r1 = norm(x - (b*bv'));
    
    % accept & repeat
    if r1 < r0
        U = b;
        V = bv;
        fprintf('Iteration: %d | Improvement df = %d\n',i,(r1-r0)./r0);
    else
        fprintf('Iteration: %d | Stop (best f = %d)\n',i,r0);
        return;
    end

end


end

function V = singlemember(V)
% each column can only belong to 1 row
for i = 1:size(V,1)

    I = find(V(i,:));

    if length(I) > 1
        
        V0(i,:)  = V(i,:);
        [~,pnt]  = max(abs(V(i,:)));
        V(i,:)   = V(i,:)*0;
        V(i,pnt) = V0(i,pnt);

    end

end


end