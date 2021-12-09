% setup with /Users/matt/Documents/MATLAB/KRSP/predict_awake.m
doSave = true;
subplotMargins = [.15,.07]; % [vert, horz]

close all
% statsMult = [1 -1;-1 1;1 1;-1 -1];
statsMult = [0 -1;0 1;1 0;-1 0];
h = ff(1200,400);
fs = 14;
subplot_tight(1,3,1,subplotMargins);
op = 0.075;
useIds = [2,3,1,4,2];
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);

lns = [];
for ii = 1:size(overlapStats,1)
    theseStats = overlapStats(ii,:);
    for jj = 1:4
        xs = [theseStats(useIds(jj))*statsMult(useIds(jj),1),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),1)];
        ys = [theseStats(useIds(jj))*statsMult(useIds(jj),2),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),2)];
        lns(1) = plot(xs,ys,'color',[0,0,0,op]);
        hold on;
    end
end

for iSeason = 1:4
    ss = ismember(mean_doys,seasonDoys(sIds(iSeason):sIds(iSeason+1)));
    theseStats = mean(overlapStats(ss,:));
    for jj = 1:4
        xs = [theseStats(useIds(jj))*statsMult(useIds(jj),1),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),1)];
        ys = [theseStats(useIds(jj))*statsMult(useIds(jj),2),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),2)];
        lns(iSeason+1) = plot(xs,ys,'color',[colors(iSeason,:),1],'linewidth',2.5);
        hold on;
    end
end
legend(lns,{'Individual','Winter','Spring','Summer','Autumn'});
legend box off
title('Nest-QB Overlap');

xlim([-0.5 0.75]);
ylim([-0.25 1]);
xticks(sort([0,xlim]));
yticks(sort([0,ylim]));
xticklabels(abs(xticks));
yticklabels(abs(yticks));
offset = 0.05;
text(-offset,max(ylim)-offset*2,'in-QB','horizontalalignment','right','fontsize',fs);
text(0.3,-offset*2,'out-AB','horizontalalignment','left','fontsize',fs);
text(-offset,min(ylim)+offset*2,'in-AB','horizontalalignment','right','fontsize',fs);
text(-0.4,offset*2,'out-QB','horizontalalignment','left','fontsize',fs);
% grid on;
xlabel('Fraction of Day');
ylabel('Fraction of Day');
set(gca,'fontsize',16);
box on;
grid on;
%%
f = fit(sq_inNestMin'*24,sq_asleepMin'*24,'poly1');

% newer code, plots each point individually
subplot_tight(1,3,2,subplotMargins);
for ii = 1:numel(sq_inNestMin)
    nrecs = numel(find(overlapMeta.squirrelId{ii} == [overlapMeta.squirrelId{:}]));
    x = sq_inNestMin(ii)*24;
    y = sq_asleepMin(ii)*24;
    iSeason = overlapMeta.meanSeason{ii};
    useMarkerSize = 7;
    if nrecs == 1
        useMarker = '.';
        useMarkerSize = 15;
    elseif nrecs == 2
        useMarker = '^';
    else
        useMarker = 'o';
    end
    plot(x,y,useMarker,'color',colors(iSeason,:),'markersize',useMarkerSize,'linewidth',1);
    hold on;
end

% older code, plots points for each season
season_residuals = {};
% code for interaction
for iSeason = 1:4
    ss = ismember(mean_doys,seasonDoys(sIds(iSeason):sIds(iSeason+1)));
    x = sq_inNestMin(ss)*24;
    y = sq_asleepMin(ss)*24;
    season_residuals{iSeason} = f(x)'-y;
%     plot(x,y,'.','color',colors(iSeason,:),'markersize',15);
%     hold on;
end

xlim([6 24]);
xticks(6:2:24);
ylim([6 26]);
yticks(6:2:24);

% f(x) = p1*x + p2
lns = [];
lns(1) = plot(NaN,NaN,'k.','markersize',15);
lns(2) = plot(xlim,[f.p1*min(xlim)+f.p2,f.p1*max(xlim)+f.p2],'k-');
ci = confint(f);
lns(3) = plot(xlim,[ci(1,1)*min(xlim)+ci(1,2),ci(1,1)*max(xlim)+ci(1,2)],'k:');
plot(xlim,[ci(2,1)*min(xlim)+ci(2,2),ci(2,1)*max(xlim)+ci(2,2)],'k:');

[r,p] = corr(sq_inNestMin'*24,sq_asleepMin'*24);
title('QB vs. In Nest');
xlabel('In Nest (hrs/day)');
ylabel('QB (hrs/day)');
grid on;
set(gca,'fontsize',16);
legend(lns,{'Squirrel Mean','Linear Fit','95% Confidence'},'location','northwest');
legend box off;

text(mean(xlim),min(ylim)+2,sprintf('r = %1.2f, p = %1.2e',r,p),'horizontalalignment','center','fontsize',14);
text(mean(xlim),min(ylim)+1,sprintf('QB = %1.2f × in nest + %1.2f',f.p1,f.p2),'horizontalalignment','center','fontsize',14);

%% season_residuals
resBins = linspace(-4,4,20);
subplot_tight(1,3,3,subplotMargins);
for iSeason = 1:4
    counts = histcounts(season_residuals{iSeason},resBins) ./ numel(season_residuals{iSeason});
    x = linspace(min(resBins),max(resBins),numel(counts));
    hold on;
    xline(sum(x.*counts),'color',colors(iSeason,:),'linewidth',2);
    plot(x,counts,'-','linewidth',4,'color',colors(iSeason,:));
end
hold off;
xlim([min(resBins),max(resBins)]);
ylabel('Fraction of Squirrels');
grid on;
set(gca,'fontsize',16);
title('QB vs. In Nest Residuals');
xlabel('Residual (hrs)');
ylim([0 0.35]);
text(min(xlim)+0.05,0.335,'\leftarrow more QB in nest','horizontalalignment','left',...
    'verticalalignment','middle','fontsize',14);
text(max(xlim)-0.05,0.315,'less QB in nest \rightarrow ','horizontalalignment','right',...
    'verticalalignment','middle','fontsize',14);

% addFigureLabels(h);
% setFig('','',2); % not sure if this is needed?
set(h,'PaperPositionMode','auto');

if doSave
%     print(gcf,'-painters','-depsc',fullfile(exportPath,'nestSleepOverlap.eps'));
    saveas(h,fullfile(exportPath,'nestSleepOverlap.eps'),'epsc');
    saveas(h,fullfile(exportPath,'nestSleepOverlap.jpg'),'jpg');
    close(h);
end