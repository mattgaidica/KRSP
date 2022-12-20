if do
    T = readtable('/Volumes/GAIDICASSD/KRSP/KRSP Axy Data/Emily Studds Axy/April 2015 - Week 1/B10_Apr23_2015.csv');
    do = 0;
end
savePath = '/Users/matt/Dropbox (Personal)/Presentations/2022 KRSP';

odba = T.odba;
[temp,unshiftedNest] = getTempAndNest(T.temp,120); % perform k-means on temp, remove repeating vals
filteredNest = filterNest(unshiftedNest,odba,temp); % remove unprobable transitions
% rmShortNest = removeShortTransitions(filteredNest,10); % optional
% nest = fixTempDelay(rmShortNest,odba,temp); % re-align nest

%% find all nest trans, plot perievent temp/nest
perieventNest = removeShortTransitions(nest,60*10); % optional
nMin = 60;
windowHalfSize = 60*nMin;
diff_nest = diff(perieventNest);
enterLocs = find(diff_nest==1);
exitLocs = find(diff_nest==-1);
enterNestTemp = NaN(1,windowHalfSize*2);
exitNestTemp = NaN(1,windowHalfSize*2);

enterNest = NaN(numel(enterLocs),windowHalfSize*2);
exitNest = NaN(numel(exitLocs),windowHalfSize*2);
jj = 0;
for ii = 1:numel(enterLocs)
    useRange = enterLocs(ii)-windowHalfSize:enterLocs(ii)+windowHalfSize-1;
    if useRange(1) > 0 && useRange(end) < numel(temp)
        jj = jj + 1;
        enterNestTemp(jj,:) = temp(useRange);
        enterNest(jj,:) = perieventNest(useRange);
    end
end
jj = 0;
for ii = 1:numel(exitLocs)
    useRange = exitLocs(ii)-windowHalfSize:exitLocs(ii)+windowHalfSize-1;
    if useRange(1) > 0 && useRange(end) < numel(temp)
        jj = jj + 1;
        exitNestTemp(jj,:) = temp(useRange);
        exitNest(jj,:) = perieventNest(useRange);
    end
end
t = linspace(-nMin,nMin,size(enterNestTemp,2));
% close all;
ff(1200,400);
titleLabels = {'Enter Nest','Exit Nest'};
tempData = {enterNestTemp,exitNestTemp};
nestData = {enterNest,exitNest};
for ii = 1:2
    subplot(1,2,ii);
    plot(t,tempData{ii}');
    hold on;
    ln1 = plot(t,mean(tempData{ii}),'color',[1 0 0 0.4],'linewidth',4);
    if ii == 1
        useylim = ylim;
    else
        ylim(useylim);
    end
    ylabel('Temp (C)');

    yyaxis right;
    plot(t,nestData{ii}','-','color',[0 0 0 0.4],'linewidth',2);
    yticks([0 1]);
    yticklabels({'Out of Nest','In Nest'});
    set(gca,'ycolor','k');
    set(gca,'fontsize',14);
    ylim([-1 2]);
    title(sprintf("%s (n=%i)",titleLabels{ii},size(nestData{ii},1)));
    xlabel('Time (min)');
    grid on;
    legend(ln1,'avg temp');
end
saveas(gcf,fullfile(savePath,'enter-exit-nest-temp-data.jpg'));
%% ^plot individually
ff(1000,400);
for ii = 1:size(exitNestTemp,1)
    plot(t,exitNestTemp(ii,:));
    ylim([20 33]);
    grid on;
    hold on;
    yline(27.7,'r');
end

%% temp histogram
colors = lines(5);
% close all
ff(1200,400);
subplot(121);
histogram(temp,100,'facecolor','k','edgecolor','none');
grid on;
set(gca,'fontsize',14);
ylabel('Samples');
xlabel('Temp (C)');
title('Temp Histogram');

[IDX,C] = kmeans(temp,2);
subplot(122);
histogram(temp(nest==1),50,'edgecolor','none','facecolor','k','facealpha',0.75);
hold on;
histogram(temp(nest==0),50,'edgecolor','none','facecolor',colors(5,:),'facealpha',0.75);
grid on;
set(gca,'fontsize',14);
ylabel('Samples');
xlabel('Temp (C)');
title('In/Out Nest by K-means');
legend({'In Nest','Out of Nest'},'location','northwest');
saveas(gcf,fullfile(savePath,'temp-histograms-k-means.jpg'));

%% plot axy/temp/nest subsection
startSample = 1+60*60* 48;
useSamples = 60*60*48; % hours
useRange = startSample:startSample+useSamples-1;
tempRange = temp(useRange);
oldNestRange = strcmp(T.Nest(useRange),'Nest');
odbaRange = odba(useRange);
nestRange = nest(useRange);
% nestRange = strcmp(T.Nest(useRange),'Nest');

colors = lines(5);
t = linspace(0,useSamples/60,useSamples); % minutes
% close all
ff(1200,400);
ylabel('raw axy');
plot(t,odbaRange,'k','linewidth',2);
ylabel('ODBA');
set(gca,'fontsize',14);

yyaxis right;
plot(t,nestRange,'color',colors(5,:),'linewidth',3);
hold on;
% plot(t,oldNestRange,'-','color',[0.2 1 0.5 0.5],'linewidth',2);
plot(t,normalize(tempRange,'range'),'r-','linewidth',2);
set(gca,'ycolor',colors(5,:));
set(gca,'fontsize',14);
ylim([-.1 1.1]);
yticks([0 1]);
yticklabels({'Out of Nest','In Nest'});
xlim([min(t),max(t)]);
xlabel('Time (min)');
grid on;
legend({'ODBA','New Nest','Temp'})
title('December 2015');
% saveas(gcf,fullfile(savePath,'axy-temp-data-overview.jpg'));
