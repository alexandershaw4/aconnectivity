function [innet,netout,clustersize,hub] = definesubnet(net,v,ET,pc)
% from a functional connectivity network, this routine will identify the
% node with highest degree (most connected), then the whole network
% attached to it s.t. assumption at least 'ET' number of edges must touch
% for a node to be in the network.
%
% usage: [innet,netout] = definesubnet(net,[v],[ET],[percent])
%
% where: net = n*n FC matrix and other inputs are optional:
%        v = (optional) coordinates of nodes in network (uses this to
%        weight the FC matrix during thresholding).
%        ET = min number of edges must touch to be network (default = 2)
%        percent = %-variance to retain when thresholding matrix [*]
%
% [*] this runs a small iterative optimismation to find the top abs max n-%
% threshold that coresponds to retaining percent% variance.
%
% AS2023

if nargin < 3 || isempty(ET)
    ET = 2;
end
if nargin < 4 || isempty(pc)
    pc = 0.8;
end

% preprep: symmetric and no diagonal
net  = net .* ~eye(length(net));
net  = (net + net') / 2;
netx = net;
%wt   = aconnectivity.QtoGauss(ones(length(net),1),8);
%net  = net.*wt;

if pc < 1
    if nargin > 1 && ~isempty(v)
        % use distances between physical nodes
        D = aconnectivity.cdist(v,v);
        ND = 1-(D./max(D(:)));
        n = aconnectivity.threshfind(net.*ND,pc)./ND;
    else
        % or just use functinoal matrix
        n = aconnectivity.threshfind(net,pc);
    end
else
    fprintf('Skipping thresholding\n');
    n = net;
end

n = (n + n')./2;
n = denan(n);

% hub index
i   = aconnectivity.hubi(n);
hub = i;

% fixed to not only be nextdoor
[innet,st] = aconnectivity.identify(n,i,ET,0);

% what to do if empty i.e. hub wrong component
if isempty(innet)
    search = true; cnt = 0;
    while search
        n(i,:) = 0;
        n(:,i) = 0;
    
        i = aconnectivity.hubi(n);
    
        [innet,st] = aconnectivity.identify(n,i,ET);

        cnt = cnt + 1;

        if cnt >= 1000
            search = false;
        end
    end

end

innet = [i st; innet];

iterate = 1; cnt = 0;
while iterate && cnt < 1000
    cnt = cnt + 1;
    l0 = size(innet,1);
    innet = aconnectivity.identify(n,innet,ET);
    l1 = size(innet,1);

    fprintf('it: %d clustersize = %d (%d %%-change)\n',cnt,l1,100*(l1-l0)./l0);
    clustersize= l1;

    if l1 == l0;
        iterate = 0;
    end
end

netout = net*0;
for i = 1:length(innet)
    netout(innet(i,1),innet(i,2)) = n(innet(i,1),innet(i,2));
end

netout = (triu(netout) + triu(netout')) + (tril(netout) + tril(netout)') ;
netout = netout./2;

end

% function [innet,ind] = identify(n,i,T)    
% 
%     innet = [];
% 
%     if nargin < 3 || isempty(T)
%         T = 1;
%     end
% 
%     for j = 1:length(i)
% 
%         if size(i,2) == 1
%             % seed
%             [~,ind] = max(abs(n(i(j,1),:)));
%         else
%             ind = i(j,2);
%         end
% 
%         % neighbours to seed
%         l = ineighb([i(j) ind]);
% 
%         % check neighbours exist e.g. at edge of grid
%         l = checkneighb(l,n);
% 
%         % check is nonzero (is connected)
%         %ic = find(diag(n(l(:,1),l(:,2))));
% 
%         for ik = 1:size(l,1)
%             ic(ik) = n(l(ik,1),l(ik,2));
%         end
%         ic = find(ic);
% 
% 
%         if length(ic) >= T
%             innet = [innet; l(ic,:)];
%         end
% 
%     end
% 
%     innet = unique(innet, 'rows');
% 
% end
% function n = ineighb(x)
%     % adjacent neighbours on a square matrix
%     n = [ x - [1 0];
%           x - [1 1];
%           x - [0 1];
%           x + [1 0];
%           x + [1 1];
%           x + [0 1];
%           x + [-1 1];
%           x + [1 -1] ];
% 
% 
% end

% function l = checkneighb(l,n);
%     % check coords exist and dont stray outside edges of matrix
%     m = size(n,1);
%     g = 1:size(l,1);
% 
%     for i = 1:size(l,1)
%         if l(i,1) > m || l(i,2) > m
%             g(i) = 0;
%         end
%     end
%     g = g(g>0);
%     l = l(g,:);
% 
% end

% function n = threshfind(net,p)
% % simple gradient descent to find threshold value to retain n-% variance of
% % input matrix
% 
% k = 0.8;
% n = thresh(net,k);
% r = corr(net(:),n(:)).^2;
% e = (r-p)./p;
% 
% iterate = 1;
% it = 0;
% while iterate
%     it = it + 1;
%     if r > p;
%         k = k - 0.1*e;
%     elseif r < p
%         k = k + 0.1*e;
%     end
% 
%     n = thresh(net,k);
%     r = corr(net(:),n(:)).^2;
% 
%     fprintf('it: %d: r^2 = %d (e = %d)\n',it,r,e);
% 
%     if r < (p*1.05) && r > (p*0.95)
%         iterate = 0;
%     end
% 
% end
% 
% end


% function i = hubi(net)
% 
% net = net.*~eye(length(net));
% [~,i] = max(sum(abs(net)));
% 
% end

% function net = thresh(net,np)
% 
% net = net.*~eye(length(net));
% n = net(:);
% 
% [~,I] = maxpoints((abs(n)),np*(length(find(n)))); % get indices
% b = find(~ismember(1:(length(n)),I));
% n(b) = 0;
% net = reshape(n,[length(net),length(net)]);
% 
% end