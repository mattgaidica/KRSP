close all
ff(1200,400);
nestTitles = {'Just Exited Nest','Just Entered Nest'};
for iState = 1:2
    Z = [];
    Zdoys = [];
    ii = 0;
    for iSq = 1:numel(unSqRow)
        useRows = find(nestMeta(:,1) == stateKey(iState) & nestMeta(:,3) == unSqRow(iSq));
        if ~isempty(useRows)
            ii = ii + 1;
            theseOdba = nestOdba(useRows,:);
            meanOdba = nanmean(theseOdba);
            Z(ii,:) = meanOdba;
            Zdoys(ii) = nestMeta(useRows(1),2);
        end
    end
    subplot(1,2,iState);
    [sortedDoys,k] = sort(Zdoys);
    s = surf(1:numel(meanOdba),1:ii,Z(k,:),'FaceAlpha',1);
    s.EdgeColor = 'none';
    colormap(magma);
    title(nestTitles{iState});
    xlabel('Time (min)');
    ylabel('Squirrel/DOY');
    zlabel('ODBA');
    set(gca,'fontsize',12);
    caxis([0.075 1])
    ylim([1 iSq]);
    xlim([1 41]);
end
view(46,37);
saveas(gcf,'/Users/matt/Desktop/entry-exit-surf-plot.jpg');

for iState = 1:2
    subplot(1,2,iState);
    view(90,-90);
end
saveas(gcf,'/Users/matt/Desktop/entry-exit-surf-plot-flat.jpg');

%%
ff(400,300);
plot(sortedDoys,'k','linewidth',2);
xlim([1 numel(sortedDoys)]);
ylim([-3 max(sortedDoys)]);
xlabel('Squirrel');
ylabel('DOY');
hold on;
for iSeason = 1:4
    plot(useDoys{iSeason},-2,'.','color',seasonColors(iSeason,:),'markerSize',15);
end
set(gca,'fontsize',14);
saveas(gcf,'/Users/matt/Desktop/doy-vs-squirrel.jpg');

%%
close all
ff(1200,400);
nestTitles = {'Just Exited Nest','Just Entered Nest'};
for iState = 1:2
    Z = [];
    Zdoys = [];
    ii = 0;
    for iSq = 1:numel(unSqRow)
        useRows = find(nestMeta(:,1) == stateKey(iState) & nestMeta(:,3) == unSqRow(iSq));
        if ~isempty(useRows)
            ii = ii + 1;
            theseOdba = nestOdba(useRows,:);
            meanOdba = nanmean(theseOdba);
            Z(ii,:) = gradient(meanOdba);
            Zdoys(ii) = nestMeta(useRows(1),2);
        end
    end
    subplot(1,2,iState);
    [sortedDoys,k] = sort(Zdoys);
    s = surf(1:numel(meanOdba),1:ii,Z(k,:),'FaceAlpha',1);
    s.EdgeColor = 'none';
    colormap(magma);
    title(nestTitles{iState});
    xlabel('Time (min)');
    ylabel('Squirrel/DOY');
    zlabel('ODBA');
    set(gca,'fontsize',12);
    caxis([-.1 .1]);
    ylim([1 iSq]);
    xlim([1 41]);
    view(90,-90);
    c = colorbar;
    c.Label.String = 'ODBA';
end
saveas(gcf,'/Users/matt/Desktop/entry-exit-surf-gradient.jpg');