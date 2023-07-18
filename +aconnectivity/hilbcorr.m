function R = hilbcorr(x,fq,fs)

    for i = 1:size(x,1)
        fx(i,:) = atcm.fun.bandpassfilter(x(i,:),fs,[fq(1) fq(2)]);
    end
    
    hx = abs(hilbert(fx));
    R  = corr(hx');

end