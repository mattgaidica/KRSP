% generate RITable: /Users/matt/Documents/MATLAB/KRSP/xcorrAsleep.m
close all
ff(450,1000);

useIds = find(RITable.is_mast==1 & ~isnan(RITable.midden_cones_diff) & RITable.is_mast==1);
xs = {normalize(RITable.qb(useIds)),normalize(RITable.qb_day(useIds)),normalize(RITable.qb_night(useIds)),normalize(RITable.RI_odba(useIds)),...
    normalize(RITable.age(useIds))};
y = RITable.midden_cones_diff(useIds);
idsSex = RITable.sex(useIds);
colors = lines(5);
plotLabels = {'QB','QB_{day}','QB_{night}','RI','Age'};
for iPlot = 1:numel(plotLabels)
    subplot(numel(plotLabels),1,iPlot);
    [f,gof] = fit(xs{iPlot},y,'poly1');
    [r,p] = corr(xs{iPlot},y);
    scatter(xs{iPlot},y,45,colors(idsSex+4,:),'filled');
    hold on;
    plot(xs{iPlot},f(xs{iPlot}),'r-');
    text(xs{iPlot},y,compose(' %s',num2str(RITable.squirrel_id(useIds))),'color',repmat(0.6,[1,3]));
    xlabel([plotLabels{iPlot},' (Z)']);
    ylabel('NMC (Z)');
    title(['New Midden Cones vs. ',plotLabels{iPlot}]);
    grid on;
    set(gca,'fontsize',12);
    ylim([-2 2]);
    yticks([-2 0 2]);
    xticks(floor(min(xlim)):ceil(max(xlim)));
    text(max(xlim),min(ylim),sprintf('r^2_{adj}=%1.3f, p=%1.3f',gof.adjrsquare,p),'horizontalalignment','right','verticalalignment','bottom','color','k');
    hold on;
    ln1 = plot(-100,-100,'.','markersize',25,'color',colors(4,:));
    ln2 = plot(-100,-100,'.','markersize',25,'color',colors(5,:));
    legend([ln1,ln2],{'Female','Male'},'location','eastoutside','fontsize',12,'box','off');
end

% print(gcf,'-painters','-depsc',fullfile(exportPath,'NewMiddenCones_vs_X.eps')); % required for vector lines
saveas(gcf,fullfile(exportPath,'NewMiddenCones_vs_X.jpg'),'jpg');

%%
close all
ff(700,400);

useIds = find(RITable.is_mast==1 & ~isnan(RITable.midden_cones_diff) & RITable.is_mast==1);
xs = {normalize(RITable.qb_night(useIds)),normalize(RITable.qb(useIds)),normalize(RITable.qb_day(useIds)),normalize(RITable.RI_odba(useIds)),...
    normalize(RITable.age(useIds))};
y = RITable.midden_cones_diff(useIds);
idsSex = RITable.sex(useIds);
colors = lines(5);
plotLabels = {'QB_{night}','QB','QB_{day}','RI'}; % ,'Age (Z)'
plotPlaces = {[1:3,5:7,9:11],4,8,12};
rows = 3;
cols = 4;
for iPlot = 1:numel(plotLabels)
    if iPlot == 1
        ms = 50;
    else
        ms = 25;
    end
    subplot(rows,cols,plotPlaces{iPlot});
    [f,gof] = fit(xs{iPlot},y,'poly1');
    [r,p] = corr(xs{iPlot},y);
    scatter(xs{iPlot},y,ms,colors(idsSex+4,:),'filled');
    hold on;
    plot(xs{iPlot},f(xs{iPlot}),'r-');
    
    if iPlot == 1
        text(xs{iPlot},y,compose(' %s',num2str(RITable.squirrel_id(useIds))),'color',repmat(0.6,[1,3]),'fontsize',12);
        xlabel([plotLabels{iPlot},' (Z)']);
        ylabel('New Midden Cones (Z)');
        title(['New Midden Cones vs. ',plotLabels{iPlot}]);
        grid on;
        set(gca,'fontsize',14);
        pos = get(gca,'position');
        set(gca,'position',pos.*[1,1,1,0.9]+[0 0.0485 0 0]);
        hold on;
        ln1 = plot(-100,-100,'.','markersize',45,'color',colors(4,:));
        ln2 = plot(-100,-100,'.','markersize',45,'color',colors(5,:));
        legend([ln1,ln2],{sprintf('Female (n = %i)',sum(~idsSex)),sprintf('Male (n = %i)',sum(idsSex))},'location','northwest','fontsize',14,'box','off');
    else
        xlabel([plotLabels{iPlot},' (Z)']);
        grid on;
        set(gca,'fontsize',11);
    end
    ylim([-2 2]);
    yticks([-2 0 2]);
    xticks(floor(min(xlim)):ceil(max(xlim)));
    text(max(xlim),min(ylim),sprintf('r^2_{adj}=%1.3f, p=%1.3f',gof.adjrsquare,p),'horizontalalignment','right','verticalalignment','bottom','color','k');
end
xs = [-19,-3.2892];
ys = 16.5;
text(xs(1),ys,'A','fontsize',24);
text(xs(2),ys,'B','fontsize',24);

saveas(gcf,fullfile(exportPath,'NewMiddenCones_vs_X_grid.jpg'),'jpg');