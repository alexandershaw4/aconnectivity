function n = threshfind(net,p)
% simple gradient descent to find threshold value to retain n-% variance of
% input matrix

k = 0.8;
n = aconnectivity.thresh(net,k);
r = corr(net(:),n(:)).^2;
e = (r-p)./p;

iterate = 1;
it = 0;
while iterate
    it = it + 1;
    if r > p;
        k = k - 0.1*abs(e);
    elseif r < p
        k = k + 0.1*abs(e);
    end

    n = aconnectivity.thresh(net,k);
    r = corr(net(:),n(:)).^2;
    e = (r-p)./p;

    fprintf('it: %d: r^2 = %d (e = %d)\n',it,r,e);

    if r < (p*1.05) && r > (p*0.95)
        iterate = 0;
    end

end

end