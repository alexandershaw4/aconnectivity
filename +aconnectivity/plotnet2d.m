function plotnet2d(net,labels)

[xo,yo] = aconnectivity.circle(length(net),3);
[xi,yi] = aconnectivity.circle(length(net),2.94);

% generate boundaries plot
figure('Position',[1741         115        1457        1375]);
plot(xo,yo,'color',[.2 .6 1],'linewidth',4); hold on;
plot(xi,yi,'color',[.2 .6 1],'linewidth',3);
axis square;
scatter(xo,yo,30,[.4 .4 .4],'filled');
scatter(xi,yi,30,[.4 .4 .4],'filled');
axis off;

cut = round(length(net)./2);
flip = ones(1,length(net));
%flip(cut+1:end) = 180;

% add labels
for i = 1:length(labels)
    t(i) = text(xo(i),yo(i),strrep(labels{i},'_',' '), ...
        'HorizontalAlignment','left','VerticalAlignment', 'bottom','rotation',i*flip(i));
end

% add network
pos = [xi(:) yi(:)];
A = net;
node1 = []; node2 = []; strng = [];
for i = 1:length(A)
    [ix,iy,iv] = find(A(i,:));
    
    if ~isempty(ix)
        conns = max(length(ix),length(iy));
        for nc = 1:conns
            node1 = [node1; pos(i(1),:)];
            node2 = [node2; pos(iy(nc),:)];
            strng = [strng; iv(nc)];
        end
    end
end

for i = 1:length(node1)
    
    to = node1(i,:);
    from = node2(i,:);

    xline = linspace(to(1),from(1),20);
    yline = linspace(to(1),from(1),20);

    zh = prod(abs(to))./2;
    wl = hamming(20)*zh;
    wl  =rescale(wl,1,2);

    xline = xline(:).*wl(:);
    yline  =yline(:).*wl(:);

    line(xline,yline,'color',strng(i));

end