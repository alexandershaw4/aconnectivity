function subnetmask = genmask(n,indices)

K = length(indices);
for ik = 1:K
        mask = zeros(n);
        ind = indices{ik};

        for i = 1:size(ind,1)
            x = ind(i,1);
            y = ind(i,2);
            mask(x,y) = 1;
        end
        subnetmask{ik} = mask;    
end

end