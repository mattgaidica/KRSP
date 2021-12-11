function [T,sleepThresh] = detect_sleepWake3(T)
doPlot = false;
doCalc = false;
if doPlot
    doCalc = true;
end
sleepThresh = NaN;

if doCalc
    meanArr = [];
    minRef = linspace(0,86400,1440); % for each minute
    secArr = secDay(T.datetime); % save compute time
    for ii = 1:numel(minRef)-1
        idx = find(secArr > minRef(ii) & secArr <= minRef(ii+1));
        meanArr(ii) = mean(T.odba(idx));
    end

    [counts,binEdges] = histcounts(meanArr,'normalization','probability','BinMethod','sqrt');
    xbins = binEdges(1:end-1) + mean(diff(binEdges))/2;
    moreBy = 1000;
    counts = equalVectors(counts,moreBy);
    xbins = equalVectors(xbins,moreBy);
    counts = smoothdata(counts,'gaussian',numel(counts)*0.1); % by percent
    % smoothCounts = detrend(smoothCounts);
    TF = islocalmin(counts);
    firstTF = find(TF == 1,1,'first');
    if xbins(firstTF) < max(xbins) * 0.25
        sleepThresh = xbins(firstTF);
    end
else
    sleepThresh = 0.204; % found empirically using doPlot=1
end

% % % % tRange = 1:60*3;
% % % % sleepThresh = mean(meanArr(tRange)) + std(meanArr(tRange))*2;

% counts = histcounts(T.odba,nBins);
% sleepThresh = otsuthresh(counts);

awakeIdx = T.odba > sleepThresh;
asleepIdx = T.odba <= sleepThresh;
T.awake = logical(awakeIdx);
T.asleep = logical(asleepIdx);
T.odba_z = (T.odba - mean(T.odba)) / std(T.odba);

if doPlot && ~isnan(sleepThresh)
    awakeOdba = T.odba;
    asleepOdba = T.odba;
    asleepOdba(awakeIdx) = NaN;

    rows = 3;
    cols = 2;
    fs = 14;
    close all;
    ff(800,800);
    
    subplot(rows,cols,1);
    plot(meanArr,'k','linewidth',1);
    xlim([1 numel(meanArr)]);
    title('meanArr');
    set(gca,'fontsize',fs);
    ylabel('ODBA');
    xlabel('Minute of Day');
    title('Mean ODBA by Minute of Day');
    
    subplot(rows,cols,2);
    plot(xbins,counts,'k-','linewidth',2);
    hold on;
    ln1 = xline(sleepThresh,'r-','linewidth',2);
    xticks(unique(sort([0:ceil(max(xbins)),sleepThresh])));
    title('Mean ODBA Distribution with Local Minima');
    set(gca,'fontsize',fs);
    xlim([min(xbins) max(xbins)]);
    grid on;
    ylabel('Probability');
    xlabel('ODBA');
    legend(ln1,'AB-QB Threshold');
    
    colors = lines(5);
    sleepColor = [0 0 0];
    wakeColor = colors(5,:);

    subplot(rows,cols,3:4);
    plot(awakeOdba,'color',colors(5,:));
    hold on;
    plot(asleepOdba,'k-');
    yline(sleepThresh,'r-');
    xlim([1 numel(asleepOdba)]);
    set(gca,'fontsize',fs);
    title(sprintf('Session ODBA (%i minutes) by "QB" and "AB" (%s)',numel(T.odba),datestr(T.datetime(1),1)));
    xlabel('Sample (i.e. minute)');
    ylabel('ODBA');
    legend({'QB','AB','AB-QB Threshold'});
    hold off;
    
    subplot(rows,cols,5:6);
    plot(awakeOdba,'color',colors(5,:));
    hold on;
    plot(asleepOdba,'k-');
    yline(sleepThresh,'r-');
    xlim([1 min([numel(asleepOdba),720])]);
    set(gca,'fontsize',fs);
    title('Zoomed');
    xlabel('Sample (i.e. minute)');
    ylabel('ODBA');
    legend({'QB','AB','AB-QB Threshold'});
    hold off;
    ylim([0 0.5]);
end