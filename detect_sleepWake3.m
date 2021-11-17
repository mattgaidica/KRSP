function T = detect_sleepWake3(T)
doPlot = false;

meanArr = [];
medianArr = [];
minRef = linspace(0,86400,1440); % for each minute
for ii = 1:numel(minRef)-1
    idx = find(secDay(T.datetime) > minRef(ii) & secDay(T.datetime) <= minRef(ii+1));
    meanArr(ii) = mean(T.odba(idx));
    medianArr(ii) = median(T.odba(idx));
end

nBins = 24;
[counts,bins] = histcounts(medianArr,nBins);
maxBin = bins(2:end)';
idx = kmeans(counts', 2);
if mean(maxBin(idx==1)) > mean(maxBin(idx==2)) % k=2 is sleep
    kSleep = 2;
else % k=1 is sleep
    kSleep = 1;
end
sleepThresh = max(maxBin(idx==kSleep));
awakeIdx = T.odba > sleepThresh;
asleepIdx = T.odba <= sleepThresh;
T.awake = logical(awakeIdx);
T.asleep = logical(asleepIdx);
T.odba_z = (T.odba - mean(T.odba)) / std(T.odba);

if doPlot
    awakeOdba = T.odba;
    awakeOdba(asleepIdx) = NaN;
    asleepOdba = T.odba;
    asleepOdba(awakeIdx) = NaN;

    rows = 2;
    cols = 4;
    fs = 14;
    close all;
    ff(1200,800);
    subplot(rows,cols,1:2);
    plot(meanArr,'k-');
    set(gca,'fontsize',fs);
    title(sprintf('Mean ODBA vs. Minute of Day for Squirrel %i\n(%i minutes)',iSq,numel(T.odba)));
    ylabel('ODBA');
    xlabel('Minute of Day');
    xlim([0 1440]);

    subplot(rows,cols,3:4);
    plot(medianArr,'k-');
    set(gca,'fontsize',fs);
    title(sprintf('Median ODBA vs. Minute of Day for Squirrel %i\n(%i minutes)',iSq,numel(T.odba)));
    ylabel('ODBA');
    xlabel('Minute of Day');
    xlim([0 1440]);

    subplot(rows,cols,5);
    xlocs = diff(bins)/2+bins(1:end-1);
    colors = lines(5);
    lns = [];
    for ii = 1:numel(counts)
        if idx(ii) == kSleep
            useColor = [0 0 0];
            lnId = 1;
        else
            useColor = colors(5,:);
            lnId = 2;
        end
        lns(lnId) = bar(xlocs(ii),counts(ii),mean(diff(bins)),'grouped','facecolor',useColor);
        hold on;
    end
    legend({'Quiescent Behavior','Active Behavior'});
    set(gca,'fontsize',fs);
    xlim([min(bins) max(bins)]);
    xticks(linspace(0,max(bins),10));
    xticklabels(compose("%1.2f",xticks));
    xtickangle(-90);
    xlabel('ODBA');
    title({'Median ODBA Distribution','Clustered by kmeans (n=2)'});
    legend(lns,{'Quiescent Behavior','Active Behavior'});

    subplot(rows,cols,6:8);
    plot(asleepOdba,'k-');
    hold on;
    plot(awakeOdba,'-','color',colors(5,:));
    xlim([1 numel(asleepOdba)]);
    set(gca,'fontsize',fs);
    title(sprintf('All ODBA (%i minutes) by "QB" and "AB"',numel(T.odba)));
    xlabel('Sample (i.e. minute)');
    ylabel('ODBA');
    legend({'Quiescent Behavior','Active Behavior'});
end