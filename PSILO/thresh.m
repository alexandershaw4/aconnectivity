function net = thresh(net,np)

n = net(:);

[~,I] = maxpoints((abs(n)),np*(length(find(n)))); % get indices
b = find(~ismember(1:(length(n)),I));
n(b) = 0;
net = reshape(n,[length(net),length(net)]);

end