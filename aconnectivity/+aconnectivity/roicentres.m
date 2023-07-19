function r = roicentres(v,vi)

u = unique(vi);

for i = 1:length(u)

    I = find(vi == u(i));
    C = aconnectivity.spherefit(v(I,:));
    r(i,:) = C;

    if any(isnan(r(i,:))) || any(isinf(r(i,:)));
        r(i,:) = mean(v(I,:),1);
    end

end