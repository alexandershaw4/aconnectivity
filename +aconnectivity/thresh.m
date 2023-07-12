function net = thresh(net,np)

if size(net,1) == size(net,2)
    sym = 1;
    net = net.*~eye(length(net));
else sym = 0;
end

n = net(:);

[~,I] = aconnectivity.maxpoints((abs(n)),np*(length(find(n)))); % get indices
b = find(~ismember(1:(length(n)),I));
n(b) = 0;

if sym
    net = reshape(n,[length(net),length(net)]);
else
    net = reshape(n,size(net));
end

end