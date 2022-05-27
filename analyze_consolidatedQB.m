%% !RUNFIG! turn on figure; setup in predict_awake.m
close all
ff(1000,450);
doSave = 0;
doMast = [0,1]; % all: [0,1]
% alternatively, manually set do=1 then resave these
if numel(doMast) > 1
    load('allTransHist');
elseif sum(doMast) == 1
    load('allTransHist_mast');
else
    load('allTransHist_nonmast');
end

cols = 5;
rows = 2;
figs = {};
figs{1} = 1:3;
figs{2} = 6:8;
figs{3} = 4:5;
figs{4} = 9:10;

contQB = 5; % minutes
gcaFontSize = 12;
titleFontMult = 1.25;
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);


% !RUNFIG! MEAN AND PROBABILITY QB HEATMAP, set do=1 once to recompute
leftPanelPos = [1 1 1 1];
windowSize = 7;
seasonBins = round(linspace(1,366,60));
yearDoys = 1:366;
nDayBins = 144;

seasonMeanDoys = [];
for iSeason = 1:4
    seasonMeanDoys(iSeason) = useDoys{iSeason}(round(numel(useDoys{iSeason})/2));
end
allSunrise = [];
allSunset = [];
for iBin = 1:numel(seasonBins)
    allSunrise(iBin) = 0 - (secDay(Tss.noon(seasonBins(iBin))) - secDay(Tss.sunrise(seasonBins(iBin)))) ./ (3600); % if nBins=144
    allSunset(iBin) = (secDay(Tss.sunset(seasonBins(iBin))) - secDay(Tss.noon(seasonBins(iBin)))) ./ (3600);
end

clc
if do
    do = 0;
    uniqueSqs = unique(trans_is);
    allTransHist = NaN(2,numel(seasonBins),nDayBins);
    for iPlot = 1:2
        for iBin = 1:numel(seasonBins)
            theseDoys = circshift(1:366,windowSize-seasonBins(iBin));
            theseDoys = theseDoys(1:windowSize*2);
            transStarts = [];
            iDoyCount = 0;
            histDoyArr = [];
            all_durations = [];

            for iSq = uniqueSqs % must extract by recording to ensure durations are properly computed
                mastFilt = zeros(size(trans_to));
                if ismember(1,doMast)
                    mastFilt = mastFilt | ismember(trans_yr,[2014,2019]);
                end
                if ismember(0,doMast)
                    mastFilt = mastFilt | ismember(trans_yr,[2015:2018,2020]);
                end
                
                theseSleepTrans = find(trans_to==0 & ismember(trans_on,theseDoys) & trans_is==iSq & mastFilt);
                if numel(theseSleepTrans) > windowSize * 200
                    allEntries = find(trans_is==iSq);
                    % must end with awake transition to get time
                    if trans_to(allEntries(end)) == 0
                        theseSleepTrans = theseSleepTrans(1:end-1);
                    end
                    theseTransEnd = trans_at(theseSleepTrans+1); % awake
                    theseTransStart = trans_at(theseSleepTrans); % asleep
                    theseDurations = (theseTransEnd - theseTransStart) / 60; % make minutes

                    % fix transitions around next day (in seconds)
                    for ii = find(theseDurations < 0)
                        theseDurations(ii) = theseDurations(ii) + 1440;
                    end
                    all_durations = [all_durations theseDurations];
                    if iPlot == 1
                        transStarts = [transStarts theseTransStart(theseDurations > 0)];
                    else
                        transStarts = [transStarts theseTransStart(theseDurations > contQB)];
                    end
                end
            end
            [counts,binEdges,BIN] = histcounts(transStarts,linspace(0,86400,nDayBins+3),'normalization','probability'); %'BinMethod','sqrt'

            % use bin indexes to get mean sleep duration, overwrite counts
            if iPlot == 1
                for ii = 1:numel(counts)
                    counts(ii) = median(all_durations(BIN == ii));
                end
            end

            counts = counts(2:end-1); % rm edges
            binShift = round(((86400/2) - mean(secDay(Tss.noon(theseDoys)))) / (86400/nDayBins));
            allTransHist(iPlot,iBin,:) = circshift(counts,binShift);
            fprintf('%i/%i season bins\n',iBin,numel(seasonBins));
        end
    end
else
    disp('Transitions not recomputed (do=0)');
end

for iPlot = 1:2
    subplot(rows,cols,figs{iPlot});
    pos = get(gca,'Position');
    set(gca,'Position',pos.*leftPanelPos);
    
    xSolarNoon = linspace(-12,12,size(squeeze(allTransHist(iPlot,:,:)),2));
    pcolor(xSolarNoon,seasonBins,squeeze(allTransHist(iPlot,:,:)));
    shading interp;
    xticks(-12:2:12);
    useyticks = round(linspace(1,366,16));
    yticks(useyticks);
    
    colormap(gca,magma(6));
    set(gca,'ydir','normal');
    
    useyticklabels = compose("%i",useyticks);
    for iSeason = 1:4
        targetDoy = useDoys{iSeason}(round(numel(useDoys{iSeason})/2));
        yidx = closest(yticks,targetDoy);
        useyticklabels{yidx} = seasonLabels{iSeason};
    end
    yticklabels(useyticklabels);
    set(gca,'fontsize',gcaFontSize,'TitleFontSizeMultiplier',titleFontMult);
    xlabel("Rel. to Solar Noon (hours)");
    ylabel('Day of Year');
    c = colorbar;
    
    if iPlot == 1
        title('Median QB Duration');
        ylabel(c,'Minutes','fontsize',gcaFontSize);
        caxisauto(squeeze(allTransHist(iPlot,:,:)),1)
    else
        title('Consolidated QB Probability');
        ylabel(c,'Probability','fontsize',gcaFontSize);
        caxisauto(squeeze(allTransHist(iPlot,:,:)),1)
    end
    
    %     caxis(caxis*0.8) % amplify colors
    text(0,295,'No Data','verticalalignment','middle','horizontalalignment','center','fontsize',gcaFontSize);
    
    hold on;
    plot(allSunrise,seasonBins,'color',[1 1 1 0.6],'linewidth',4);
    plot(allSunset,seasonBins,'color',[1 1 1 0.6],'linewidth',4);
    hold off;
end


% xSolarNoon = 1:144[-12:12], seasonBins = doys(n=60), allSunrise/Sunset = times relative to xSolarNoon
seasonMedians = NaN(numel(seasonBins),1);
seasonProbs = NaN(numel(seasonBins),1);
seasonArr = NaN(numel(seasonBins),1);
for yy = 1:numel(seasonBins)
    sunriseId = closest(xSolarNoon,allSunrise(yy));
    sunsetId = closest(xSolarNoon,allSunset(yy));
    ssRange = [1:sunriseId sunsetId:numel(xSolarNoon)];
    seasonMedians(yy) = mean(allTransHist(1,yy,ssRange));
    seasonProbs(yy) = mean(allTransHist(2,yy,ssRange));
    if ~isnan(seasonProbs(yy))
        % find season
        for thisSeason = 1:4
            if any(ismember(useDoys{thisSeason},seasonBins(yy)))
                seasonArr(yy) = thisSeason;
            end
        end
    end
end
% % % % medianHeights = [];
% % % % probHeights = [];
% % % % for iSeason = 1:4
% % % %     medianHeights(iSeason) = nanmean(seasonMedians(seasonArr==iSeason));
% % % %     probHeights(iSeason) = nanmean(seasonProbs(seasonArr==iSeason));
% % % % end
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
CData = NaN(numel(seasonBins),3);
for ii = 1:size(CData,1)
    if ~isnan(seasonArr(ii))
        CData(ii,:) = colors(seasonArr(ii),:);
    end
end

useData = {seasonMedians,seasonProbs};
useTitles = {'Avg. Median QB Duration (minutes)','Avg. Consolidated QB Probability'};
useXLims = [20,.0175];
legendLoc = {'northeast','southeast'};
seasonOrder = [1,4,3,2];
textStrings = {'',''};
if numel(doMast) > 1
    textStrings = {{'Spring × All *P < 0.05','Others N.S.'},'All × All *P < 0.05'};
end
for iPlot = 1:2
    subplot(rows,cols,figs{iPlot+2});
    b = barh(useData{iPlot},'facecolor','flat','edgecolor','none');
    b.CData = CData;
    yticks([]);
    ylabel('Day of Year');
    title('night-only');
    set(gca,'fontsize',gcaFontSize);
    xlabel(useTitles{iPlot});
    ylim([1 numel(seasonBins)]);
    xlim([0 useXLims(iPlot)]);
    set(gca,'XGrid','on');
    hold on;
    lns = [];
    for iSeason = 1:4
        lns(iSeason) = bar(1,1,'facecolor',colors(seasonOrder(iSeason),:),'edgecolor','none');
        xVal = nanmean(useData{iPlot}(seasonArr==iSeason));
        plot([xVal xVal],[min(ylim) max(ylim)],'color',[colors(iSeason,:),0.5],'linewidth',2);
    end
    le = legend(lns,seasonLabels(seasonOrder),'box','off','location',legendLoc{iPlot});
    text(0.62,0.5,textStrings{iPlot},'Units', 'Normalized','verticalAlignment','middle');
    hold off;
end
xs = [-0.0338,-0.001];
ys = [149,66.63];
text(xs(1),ys(1),'A','fontsize',24);
text(xs(2),ys(1),'B','fontsize',24);
text(xs(1),ys(2),'C','fontsize',24);
text(xs(2),ys(2),'D','fontsize',24);

if numel(doMast) > 1
    filename = 'ConsolidatedQB_v2_histograms';
elseif sum(doMast) == 1
    filename = 'ConsolidatedQB_v2_histograms_mast';
else
    filename = 'ConsolidatedQB_v2_histograms_nonmast';
end
saveas(gcf,fullfile(exportPath,[filename,'.jpg']),'jpg');
%%
[p,tbl,stats] = anovan(seasonMedians,seasonArr,'display','off');
disp('Medians');
c_median = multcompare(stats,'display','off')

[p,tbl,stats] = anovan(seasonProbs,seasonArr,'display','off');
disp('Probs');
c_prob = multcompare(stats,'display','off')