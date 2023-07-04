function l = checkneighb(l,n);
    % check coords exist and dont stray outside edges of matrix
    m = size(n,1);
    g = 1:size(l,1);

    for i = 1:size(l,1)
        if l(i,1) > m || l(i,2) > m
            g(i) = 0;
        end
    end
    g = g(g>0);
    l = l(g,:);

end