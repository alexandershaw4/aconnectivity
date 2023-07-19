function V = computereducedoperator(M)

[~,ix] = find(sum(M));
V = sparse(1:length(ix),(ix),1,length(ix),size(M,2));


end