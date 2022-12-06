close all;
iRec = iRec + 1;
colors = lines(5);
ff(1200,500);
plot(sq_odba(iRec,:),'k');
ylabel('ODBA');

yyaxis right;
plot(sq_nest(iRec,:),'color',colors(5,:));
ylim([-1 2]);
yticks([0 1]);
yticklabels({'Out of Nest','In Nest'});
set(gca,'ycolor',colors(5,:));
set(gca,'fontsize',16);
xlim(size(sq_nest(iRec,:)));
xlabel('Time (s)');
title(sprintf('1-day of Axy & Nest Data (%s)',seasonLookup(sq_doys(iRec))));

saveas(gcf,sprintf('/Users/matt/Desktop/axy-overview-iRec%03d.jpg',iRec));

%%
close all
ff(1000,500);
plot(theseOdba');
ylabel('ODBA');
xline(21);
xlim([1 41]);
xticks([1 21 41]);
xticklabels({'-20','Nest Entry','+20'});
title('ODBA at +/-20 min Nest Entry');
set(gca,'fontsize',16);
ylim([0 10]);

yyaxis right;
lns1 = plot(mean(theseOdba),'k-','linewidth',4);
hold on;
lns2 = plot(median(theseOdba),'-','linewidth',4,'color',colors(4,:));
legend([lns1 lns2],{'Mean','Median'});
ylabel('ODBA');
saveas(gcf,sprintf('/Users/matt/Desktop/nest-entry-mean-median.jpg',iRec));

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

%%
close all
h1 = ff(900,900);
nestTitles = {'Just Exited Nest','Just Entered Nest'};
rows = 2;
cols = 2;
for iState = 1:2
    figure(h1);
    for iSeason = 1:4
        subplot(rows,cols,iState);
        meanOdba = squeeze(sqOdba(iState,sqSeason==iSeason,:));
        seasonOdba = nanmean(meanOdba);
        plot(seasonOdba,'color',seasonColors(iSeason,:),'linewidth',3);
        hold on;
        
        xline(preMinutes+1);
        xlim([1 preMinutes+postMinutes+1]);
        title(nestTitles{iState});
        grid on;
        xlabel('minutes');
        set(gca,'fontsize',14);
        ylabel('ODBA');
        
        subplot(rows,cols,iState+2);
        meanGradient = squeeze(sqGradient(iState,sqSeason==iSeason,:));
        seasonGradient = nanmean(meanGradient);
        plot(seasonGradient,'color',seasonColors(iSeason,:),'linewidth',3);
        hold on;
        
        xline(preMinutes+1);
        xlim([1 preMinutes+postMinutes+1]);
        title(nestTitles{iState});
        grid on;
        xlabel('minutes');
        set(gca,'fontsize',14);
        ylabel('ODBA Gradient');
    end
%     [p,anovatab,stats] = anova1(decayRates,group,'off');
%     c = multcompare(stats,'display','off');
end
saveas(gcf,'/Users/matt/Desktop/ODBA-and-ODBAgradient-by-season.jpg');