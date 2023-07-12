function i = hubi(net)

net = net.*~eye(length(net));
[~,i] = max(sum(abs(~~net)));

end