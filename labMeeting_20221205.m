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

%%
close all
h1 = ff(900,400);
nestTitles = {'Just Exited Nest','Just Entered Nest'};
seasonOdba = [];
rows = 1;
cols = 2;
for iState = 1:2
    figure(h1);
    subplot(rows,cols,iState);
    for iSeason = 1:4
        meanOdba = squeeze(sqOdba(iState,sqSeason==iSeason,:));
        seasonOdba(iSeason,:) = nanmean(meanOdba);
        plot(seasonOdba(iSeason,:),'color',seasonColors(iSeason,:),'linewidth',3);
        hold on;
    end
    xline(preMinutes+1);
    xlim([1 preMinutes+postMinutes+1]);
    title(nestTitles{iState});
    grid on;
    xlabel('minutes');
    set(gca,'fontsize',16);
    ylabel('ODBA');
    xticks([1 21 41]);
    xticklabels({'-20','Nest','+20'});
    if iState == 1
        legend(seasonLabels,'location','northwest');
    else
        legend(seasonLabels);
    end
end