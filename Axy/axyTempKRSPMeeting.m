if do
    T = readtable('/Volumes/GAIDICASSD/KRSP/KRSP Axy Data/Emily Studds Axy/April 2015 - Week 1/B10_Apr23_2015.csv');
    do = 0;
end
savePath = '/Users/matt/Dropbox (Personal)/Presentations/2022 KRSP';

%% find all nest trans, plot perievent temp/nest
nMin = 60;
windowHalfSize = 60*nMin;
nest = strcmp(T.Nest,'Nest');
temp = T.temp;
temp_filt = temp;
temp_filt(temp==22.5) = NaN;
temp_filt = inpaint_nans(temp_filt);
temp_smooth = smoothdata(temp_filt,'gaussian',60);
diff_nest = diff(nest);
enterLocs = find(diff_nest==1);
exitLocs = find(diff_nest==-1);
enterNestTemp = NaN(1,windowHalfSize*2);
exitNestTemp = NaN(1,windowHalfSize*2);

enterNest = NaN(numel(enterLocs),windowHalfSize*2);
exitNest = NaN(numel(exitLocs),windowHalfSize*2);
jj = 0;
for ii = 1:numel(enterLocs)
    useRange = enterLocs(ii)-windowHalfSize:enterLocs(ii)+windowHalfSize-1;
    if useRange(1) > 0 && useRange(end) < numel(temp_smooth)
        jj = jj + 1;
        enterNestTemp(jj,:) = temp_smooth(useRange);
        enterNest(jj,:) = nest(useRange);
    end
end
jj = 0;
for ii = 1:numel(exitLocs)
    useRange = exitLocs(ii)-windowHalfSize:exitLocs(ii)+windowHalfSize-1;
    if useRange(1) > 0 && useRange(end) < numel(temp_smooth)
        jj = jj + 1;
        exitNestTemp(jj,:) = temp_smooth(useRange);
        exitNest(jj,:) = nest(useRange);
    end
end
t = linspace(-nMin,nMin,size(enterNestTemp,2));
close all;
ff(1200,400);
titleLabels = {'Enter Nest','Exit Nest'};
tempData = {enterNestTemp,exitNestTemp};
nestData = {enterNest,exitNest};
for ii = 1:2
    subplot(1,2,ii);
    plot(t,tempData{ii}');
    ylim([20 33]);
    ylabel('Temp (C)');

    yyaxis right;
    plot(t,nestData{ii}','k-');
    yticks([0 1]);
    yticklabels({'Out of Nest','In Nest'});
    set(gca,'ycolor','k');
    set(gca,'fontsize',14);
    ylim([-1 2]);
    title(titleLabels{ii});
    xlabel('Time (min)');
    grid on;
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
histogram(temp_smooth,100,'facecolor','k','edgecolor','none');
grid on;
set(gca,'fontsize',14);
ylabel('Samples');
xlabel('Temp (C)');
title('Temp Histogram');

[IDX,C] = kmeans(temp_smooth,2);
subplot(122);
histogram(temp_smooth(IDX==2),50,'edgecolor','none','facecolor','k','facealpha',1);
hold on;
histogram(temp_smooth(IDX==1),50,'edgecolor','none','facecolor',colors(5,:),'facealpha',1);
grid on;
set(gca,'fontsize',14);
ylabel('Samples');
xlabel('Temp (C)');
title('In/Out Nest by K-means');
legend({'In Nest','Out of Nest'},'location','northwest');
% saveas(gcf,fullfile(savePath,'temp-histograms-k-means.jpg'));
%% ^histogram based on nest class (not my own k-means)
nest = strcmp(T.Nest,'Nest');
temp = T.temp;
temp_filt = temp;
temp_filt(temp==22.5) = NaN;
temp_filt = inpaint_nans(temp_filt);
temp_smooth = smoothdata(temp_filt,'gaussian',60);
ff(400,400);
histogram(temp_smooth(nest==1),50);
hold on;
histogram(temp_smooth(nest==0),50);

%% plot axy/temp/nest subsection
colors = lines(5);
startSample = 60*60*48;
useSamples = 60*60*24; % hours
useRange = startSample:startSample+useSamples-1;
tempFilt = filterTemp(T,120);
temp = tempFilt(useRange);
odba = T.odba(useRange);
% nest = strcmp(T.Nest(startSample:startSample+useSamples-1),'Nest');
% [tempDelay,nestFixed] = findTempDelay(T);
nest = nestFixed(useRange);

colors = lines(5);
t = linspace(0,useSamples/60,numel(x)); % minutes
close all
ff(1200,400);
% ylabel('raw axy');
plot(t,odba,'k','linewidth',2);
ylabel('ODBA');
set(gca,'fontsize',14);
% hold on;
% plot(t,detrend((x)));
% hold on;
% plot(t,detrend((y)));
% plot(t,detrend((z)));

yyaxis right;
plot(t,nest,'color',colors(5,:),'linewidth',2);
hold on;
plot(t,normalize(temp,'range'),'r-');
set(gca,'ycolor',colors(5,:));
set(gca,'fontsize',14);
ylim([-.1 1.1]);
yticks([0 1]);
yticklabels({'Out of Nest','In Nest'});
xlim([min(t),max(t)]);
xlabel('Time (min)');
grid on;
legend({'ODBA','Nest Class','Temp'})
title('April 2015');
% saveas(gcf,fullfile(savePath,'axy-temp-data-overview.jpg'));
