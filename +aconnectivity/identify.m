function [innet,ind] = identify(n,i,T)    
    
    innet = [];

    if nargin < 3 || isempty(T)
        T = 1;
    end

    for j = 1:size(i,1)
        
        [~,ind] = max(abs(n(i(j,1),:)));

        these = find(n(i(j),:));

        list = [repmat(i(j),[length(these),1]) these(:)];
        innet = [innet; list];

        % if size(i,2) == 1
        %     % seed
        %     [~,ind] = max(abs(n(i(j,1),:)));
        % else
        %     ind = i(j,2);
        % end
        % 
        % % neighbours to seed
        % l = aconnectivity.ineighb([i(j) ind]);
        % 
        % % or switch to all connected:
        % %rl = find(n(i(j),:));
        % %l = [rl(:) ones(length(rl),1)*i(j)];
        % 
        % 
        % % check neighbours exist e.g. at edge of grid
        % l = aconnectivity.checkneighb(l,n);
        % 
        % % check is nonzero (is connected)
        % %ic = find(diag(n(l(:,1),l(:,2))));
        % clear ic
        % 
        % for ik = 1:size(l,1)
        %     ic(ik) = n(l(ik,1),l(ik,2));
        % end
        % ic = find(ic);
        % 
        % 
        % if length(ic) >= T
        %     innet = [innet; l(ic,:)];
        % end

    end

    innet = unique(innet, 'rows');

end