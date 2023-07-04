function net = thresh(net,np)

net = net.*~eye(length(net));
n = net(:);

[~,I] = aconnectivity.maxpoints((abs(n)),np*(length(find(n)))); % get indices
b = find(~ismember(1:(length(n)),I));
n(b) = 0;
net = reshape(n,[length(net),length(net)]);

end