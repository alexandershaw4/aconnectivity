function achordplot(dataMat,labels)

labels = strrep(labels,'_',' ');

rowName = labels;
colName = labels;

dataMat = (dataMat + dataMat')./2;
i = find(sum(dataMat));


% % convert to to/from/val list
% [node1,node2,vals] = conmat2nodes(dataMat(i,i),'test','sourcemodel',[(1:length(i))' (1:length(i))']);
% 
% list = [node1(:,1) node2(:,1)];


CC = chordChart(dataMat(i,i),'rowName',rowName(i),'colName',colName(i));
CC = CC.draw();

%CC.setChordColorByMap([0,0,.8;.8,0,0])


CC.setFont('FontSize',10,'FontName','Cambria')
CC.tickState('off')

% 以下代码用来旋转标签
% The following code is used to rotate the label
textHdl=findobj(gca,'Type','Text');
for i=1:length(textHdl)
    if textHdl(i).Rotation<-90
        textHdl(i).Rotation=textHdl(i).Rotation+180;
    end
    switch true
        case textHdl(i).Rotation<0&&textHdl(i).Position(2)>0
            textHdl(i).Rotation=textHdl(i).Rotation+90;
            textHdl(i).HorizontalAlignment='left';
        case textHdl(i).Rotation>0&&textHdl(i).Position(2)>0
            textHdl(i).Rotation=textHdl(i).Rotation-90;
            textHdl(i).HorizontalAlignment='right';
        case textHdl(i).Rotation<0&&textHdl(i).Position(2)<0
            textHdl(i).Rotation=textHdl(i).Rotation+90;
            textHdl(i).HorizontalAlignment='right';
        case textHdl(i).Rotation>0&&textHdl(i).Position(2)<0
            textHdl(i).Rotation=textHdl(i).Rotation-90;
            textHdl(i).HorizontalAlignment='left';
    end
end

