function con = find_between_cluster_connections(net,subnets)

k = length(subnets);

for i = 1:k
    for j = 1:k
        if i~=j 
            ix = find(sum(subnets{i}));
            iy = find(sum(subnets{j}));

            n = 0;
            for l = 1:length(ix)
                for m = 1:length(iy)
                    n = n + 1;
                    x = ix(l);
                    y = iy(m);

                    if any(net(x,y))
                        con{i,j}(n,:) = [x y];
                        con{j,i}(n,:) = [x y];
                    end
                end
            end

        end
    end
end



end