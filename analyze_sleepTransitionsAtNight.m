%% !RUNFIG! turn on figure; setup in predict_awake.m
close all
ff(500,900);
doSave = 0;
doMast = 0; % all: [0,1]

cols = 9;
rows = 12;
figs = {};
figs{1} = 1:36;
figs{2} = 37:72;
figs{3} = 73:99;
figs{4} = 100:108;

contQB = 5; % minutes
gcaFontSize = 12;
titleFontMult = 1.25;
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);

%% !RUNFIG! (1/3) MEAN AND PROBABILITY QB HEATMAP, set do=1 once to recompute
leftPanelPos = [1.2 1 0.97 0.8];
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
seasonTicks = {'Winter','Spring','Summer','Autumn'};

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
                
                theseSleepTrans = find(trans_to==0 & ismember(trans_on,theseDoys) & trans_is==iSq & mastFilt); %  & ~ismember(trans_yr,[2014,2019])
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
%         useyticklabels{yidx} = convertStringsToChars("\textbf{" + seasonTicks{iSeason} + "}");
        useyticklabels{yidx} = seasonTicks{iSeason};
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

%% !RUNFIG! (2/3) Probability Density histograms
nBins_sd = 9; % mins
showMinutes = 60;
ssBeforeAfterPad = 30; % minutes
maxY = 0.25;
leftPanelPos = [1 1 1 0.7];

% % % % close all
% % % % h = ff(450,250);
subplot(rows,cols,figs{3});
cla(gca);
lns = [];
seasonLabels = {'Winter','Spring','Summer','Autumn'};
sleepDurations_season = {};
countsMat = [];

for iSeason = 1:4
    uniqueSqs = unique(trans_is);
    sleepDurations = [];
    transStarts = [];
    for iSq = uniqueSqs
%         theseSleepTrans = find(trans_is==iSq & trans_to==0 & ismember(trans_on,useDoys{iSeason}));
        
        mastFilt = zeros(size(trans_to));
        if ismember(1,doMast)
            mastFilt = mastFilt | ismember(trans_yr,[2014,2019]);
        end
        if ismember(0,doMast)
            mastFilt = mastFilt | ismember(trans_yr,[2015:2018,2020]);
        end
        theseSleepTrans = find(trans_is==iSq & trans_to==0 & ismember(trans_on,useDoys{iSeason}) & mastFilt);
        
        % START light/dark transitions
        if ~isempty(theseSleepTrans)
            theseDoys = trans_on(theseSleepTrans);
            meanSunrise = mean(secDay(Tss.sunrise(theseDoys))) - ssBeforeAfterPad*60;
            meanSunset = mean(secDay(Tss.sunset(theseDoys))) + ssBeforeAfterPad*60;
            temp = [];
            for ii = theseSleepTrans
                if trans_at(ii) > meanSunset || trans_at(ii) < meanSunrise
                    temp = [temp ii];
                end
            end
            if isempty(temp)
                continue;
            end
            theseSleepTrans = temp;
        else
            continue;
        end
        % END light/dark code
        
        if ~isempty(theseSleepTrans)
            allEntries = find(trans_is==iSq);
            % must end with awake transition to get time
            if trans_to(allEntries(end)) == 0
                theseSleepTrans = theseSleepTrans(1:end-1);
            end
            theseTransEnd = trans_at(theseSleepTrans+1); % awake
            % make minutes, that's the highest resolution from
            % predict_awake anyways
            theseTransStart = trans_at(theseSleepTrans); % asleep
            theseDurations = (theseTransEnd - theseTransStart) / 60;
            
            % fix transitions around next day (in seconds)
            for ii = find(theseDurations < 0)
                theseDurations(ii) = theseDurations(ii) + 1440;
            end
            sleepDurations = [sleepDurations theseDurations];
        end
    end
    
    if ~isempty(sleepDurations)
        sleepDurations_season{iSeason} = sleepDurations;
        binEdges_sd = logspace(0,2,nBins_sd+1);
        if iSeason == 1 % do once
            patch([contQB,max(binEdges_sd),max(binEdges_sd),contQB],[0,0,1,1],'black','FaceAlpha',.1,'EdgeColor','none');
            hold on;
            xline(contQB,'k:');
            text(contQB,0.02,strcat(sprintf('%i minutes',contQB),'\rightarrow'),'horizontalalignment','right','fontsize',gcaFontSize);
            text(contQB,maxY-0.02,'\leftarrow fragmented ','horizontalalignment','right','fontsize',gcaFontSize);
            text(contQB,maxY-0.02,' consolidated \rightarrow','horizontalalignment','left','fontsize',gcaFontSize);
        end
        histogram(sleepDurations,binEdges_sd,'Normalization','probability','EdgeColor',colors(iSeason,:),'DisplayStyle','Stairs','lineWidth',4,'EdgeAlpha',0.75);
        hold on;
        lns(iSeason) = plot(-1,-1,'-','lineWidth',3,'color',colors(iSeason,:));
        counts = histcounts(sleepDurations,binEdges_sd,'Normalization','probability');
        countsMat(iSeason,:) = counts;
        
        % find center of 'sleep' peak
        % % % %         overSample = 10000;
        % % % %         countsSmooth = smoothdata(equalVectors(counts,overSample),'gaussian',overSample/10);
        % % % %         binsSmooth = equalVectors(binEdges,overSample);
        % % % %         locs = peakseek(countsSmooth);
        % % % %         maxBin = binsSmooth(locs(end));
        % % % %         text(17,maxY-iSeason*0.015,strcat(sprintf('%s: %1.2f',seasonLabels{iSeason}(1:2),maxBin),'\rightarrow'),...
        % % % %         'horizontalalignment','right','fontsize',12,'color',colors(iSeason,:));
        
        %         xline(median(sleepDurations),':','lineWidth',4,'color',colors(iSeason,:));
        % % % %         xline(maxBin,':','lineWidth',2,'color',colors(iSeason,:));
        set(gca,'xscale','log');
        drawnow;
    end
end

set(gca,'fontsize',gcaFontSize,'TitleFontSizeMultiplier',titleFontMult);
pos = get(gca,'Position');
set(gca,'Position',pos.*leftPanelPos);
xlabel('QB Duration (minutes)');
ylabel('Probability');
ylim([0 maxY]);
title('QB Duration at Night');
grid on;
legend(lns,seasonLabels);
legend box off;
hold off;

% % % % if doSave
% % % %     print(gcf,'-painters','-depsc',fullfile(exportPath,'DarkQBDuration.eps')); % required for vector lines
% % % %     saveas(gcf,fullfile(exportPath,'DarkQBDuration.jpg'),'jpg');
% % % %     close(h);
% % % % end

%% (3/3) setup pMat_sd
nSurr = 100;
if ~exist('pMat_sd','var')
    clc
    pMat_sd = NaN(4,4,size(countsMat,2));
    surrDiff = [];
    for iSeason = 1:4
        for kSeason = iSeason:4
            fprintf('Seasons %i-%i\n',iSeason,kSeason);
            actualDiff = diff(countsMat([iSeason,kSeason],:));
            y = [sleepDurations_season{iSeason} sleepDurations_season{kSeason}];
            group = [zeros(size(sleepDurations_season{iSeason})) ones(size(sleepDurations_season{kSeason}))];
            for iSurr = 1:nSurr
                group = group(randperm(length(group)));
                counts_i = histcounts(y(group==0),binEdges_sd,'Normalization','probability');
                counts_k = histcounts(y(group==1),binEdges_sd,'Normalization','probability');
                surrDiff(iSurr,:) = counts_i - counts_k;
            end
            for iBin = 1:size(pMat_sd,3)
                % control for direction
                if actualDiff(iBin) > 0
                    pMat_sd(iSeason,kSeason,iBin) = sum(surrDiff(:,iBin) > actualDiff(iBin)) / nSurr;
                else
                    pMat_sd(iSeason,kSeason,iBin) = sum(surrDiff(:,iBin) < actualDiff(iBin)) / nSurr;
                end
            end
        end
    end
    disp('done');
end
%% !RUNFIG! (3/3) plot pMat
cmap = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/goblue.png',200);
% cmap = hot(1000);
cmap = [cmap;0,0,0];
% % % % close all;
% % % % h = ff(1000,100);
seasonAbbr = {'Wi','Sp','Su','Au'};
for iBin = 1:size(pMat_sd,3)
    subplot(rows,cols,figs{4}(iBin));
    cla(gca);
    ps = squeeze(pMat_sd(:,:,iBin));
    imagesc(ps,'AlphaData',~isnan(ps));
    caxis([0,0.0499]);
    colormap(gca,cmap); % flip(magma)
    title(sprintf('%1.0f-%1.0f',binEdges_sd(iBin),binEdges_sd(iBin+1)));
    xticks(1:4);
    xtickangle(-90);
    yticks(1:4);
    xticklabels(seasonAbbr);
    yticklabels(seasonAbbr);
    set(gca,'fontsize',gcaFontSize-4);
    set(gca,'ydir','normal')
    pos = get(gca,'Position');
    set(gca,'Position',pos.*[1 1 0.9 0.6]+[0 -0.03 0 0]);
end
cb = cbAside(gca,'p','k');
cb.TickLabels = [0 0.05];
cb.FontSize = gcaFontSize-2;

if doSave
    %         print(gcf,'-painters','-depsc',fullfile(exportPath,'ConsolidatedQB.jpg')); % required for vector lines
    saveas(gcf,fullfile(exportPath,'ConsolidatedQB.jpg'),'jpg');
    close(gcf);
end
