% setup: /Users/matt/Documents/MATLAB/KRSP/predict_awake.m
% trans_at = [trans_at secDay(Tawake.datetime)'];
% trans_to = [trans_to Tawake.awake'];
% trans_on = [trans_on day(Tawake.datetime,'dayofyear')'];

%% SOL_T table setup (1/3), see also: /Users/matt/Documents/MATLAB/KRSP/analyze_circCorrSleep.m
close all
nSmooth = 10; % minutes
alpha = 0.2;
x = linspace(-3,3,1440);
mypdf = normalize(normpdf(x,0,1),'range');

T_SOL = table;
warning ('off','all');
iRow = 0;

doSave = 1;
doDebug = 0;
if doDebug
    debugPath = '/Users/matt/Downloads/debug';
    ms = 35;
    rows = 2;
    cols = 1;
    close all
    psSec = 0.00; % sec
    mnLabels = {'Sunrise','Sunset'};
end
 
for iRec = 1:size(sq_asleep,1)
    fprintf('iRec%i\n',iRec);
    if doDebug
        ff(1200,800);
    end
    for iSun = 1:2
        if doDebug
            subplot(rows,cols,iSun);
        end
        if iSun == 1
            shiftedAsleep = circshift(sq_asleep(iRec,:),720-round(mean(secDay(Tss.sunrise(Tss.doy == sq_doys(iRec))),1)/60));
        else
            shiftedAsleep = circshift(mean(sq_asleep(iRec,:),1),720-round(mean(secDay(Tss.sunset(Tss.doy == sq_doys(iRec))),1)/60));
        end
        
        if std(shiftedAsleep) == 0
            fprintf('skipping std rec%i\n',iRec);
            continue;
        end
        asleepNorm = normalize(imgaussfilt(shiftedAsleep,nSmooth,'Padding','circular'),'range');
        if iSun == 1
            awakeNorm = -asleepNorm + 1;
        else
            awakeNorm = asleepNorm;
        end
        locs = [];
        adjAlpha = 1-alpha;
        while isempty(locs)
            [locs,pks] = peakseek(awakeNorm.*mypdf,nSmooth,adjAlpha);
            adjAlpha = adjAlpha - 0.01;
        end
        
        if doDebug
            ln1 = plot(asleepNorm,'k','linewidth',3);
            hold on;
            ln2 = plot(mypdf,':','color',repmat(0.15,[1,3]));
            ln3 = plot(asleepNorm.*mypdf,'color',repmat(0.15,[1,3]));
            xline(720);
            yline(alpha,'r--'); yline(1-alpha,'r--');
            xlim([1 1440]);
            xticks([1,720,1440]);
            xticklabels({'-720','0','+720'});
            ylabel('QB (norm.)');
            ylim([0 1]);
            yticks([0,alpha,1-alpha,1]);
            xlabel('Time (min)');
            set(gca,'fontsize',14);
            title(sprintf('Rel. to %s, iRec = %04d',mnLabels{iSun},iRec));
            legend([ln1,ln2,ln3],{'QB','normpdf','QB × normpdf'},'Autoupdate','off');
            m1 = [];
            m2 = [];
        end
        
        if iSun == 1
            [~,awakeIdx] = closest(locs,720); % break a tie
            while asleepNorm(awakeIdx) < alpha
                if doDebug
                    delete(m1);
                    m1 = plot(awakeIdx,asleepNorm(awakeIdx),'g.','markersize',ms);
                    drawnow;
                    pause(psSec);
                end
                awakeIdx = awakeIdx - 1; % backtrack index until alpha
            end
            asleepIdx = awakeIdx;
            while asleepIdx > 1
                if doDebug
                    delete(m2);
                    m2 = plot(asleepIdx,asleepNorm(asleepIdx),'r.','markersize',ms);
                    drawnow;
                    pause(psSec);
                end
                asleepIdx = asleepIdx - 1;
                asleepAlpha = asleepNorm(asleepIdx);
                if asleepAlpha > 1-alpha
                    break;
                end
            end
            if doDebug
                plot([asleepIdx,awakeIdx],[asleepNorm(asleepIdx),asleepNorm(awakeIdx)],'color',repmat(0.5,[1,4]),'linewidth',10);
                drawnow;
            end
        else
            [~,asleepIdx] = closest(locs,720);
            while asleepNorm(asleepIdx) > 1 - alpha
                if doDebug
                    delete(m1);
                    m1 = plot(asleepIdx,asleepNorm(asleepIdx),'r.','markersize',ms);
                    drawnow;
                    pause(psSec);
                end
                asleepIdx = asleepIdx - 1;
            end
            awakeIdx = asleepIdx;
            while awakeIdx > 1
                if doDebug
                    delete(m2);
                    m2 = plot(awakeIdx,asleepNorm(awakeIdx),'g.','markersize',ms);
                    drawnow;
                    pause(psSec);
                end
                awakeIdx = awakeIdx - 1;
                asleepAlpha = asleepNorm(awakeIdx);
                if asleepAlpha < alpha
                    break;
                end
            end
            if doDebug
                plot([asleepIdx,awakeIdx],[asleepNorm(asleepIdx),asleepNorm(awakeIdx)],'color',repmat(0.5,[1,4]),'linewidth',10);
                drawnow;
            end
        end
        if doDebug
%             plot([asleepIdx,awakeIdx],[0.5,0.5],'k-','linewidth',10);
            arrow([asleepIdx,0.5],[awakeIdx,0.5],'Length',5,'Ends',iSun);
            text(max([asleepIdx,awakeIdx])+10,0.5,sprintf('SOL = %i mins',abs(asleepIdx - awakeIdx)),'HorizontalAlignment','left','fontsize',14);
        end

        iRow = iRow + 1;
        T_SOL.iRec(iRow) = iRec;
        if iSun == 1
            T_SOL.isSunrise(iRow) = 1;
        else
            T_SOL.isSunrise(iRow) = 0;
        end
        T_SOL.SOL(iRow) = abs(asleepIdx - awakeIdx);
        T_SOL.asleepIdx(iRow) = asleepIdx;
        T_SOL.awakeIdx(iRow) = awakeIdx;
        for iSeason = 1:4
            if ismember(sq_doys(iRec),useDoys{iSeason})
                 T_SOL.season(iRow) = iSeason;
            end
        end
        T_SOL.is_mast(iRow) = ismember(sq_years(iRec),[2014,2019]);
% % % %         T_SOL.asleepData(iRow) = {asleepNorm};
    end
    if doSave && doDebug
        saveas(gcf,fullfile(debugPath,sprintf('iRec%04d.jpg',iRec)));
        close gcf;
    end   
end
warning ('on','all');
writetable(T_SOL,'T_SOL');
%% load SOL_T cleaning and seasonal histograms (2/3)
T_SOL = readtable('T_SOL');
origSz = size(T_SOL,1);

T_SOL(T_SOL.awakeIdx == 1,:) = [];
T_SOL(T_SOL.asleepIdx == 1,:) = [];

fnSize = size(T_SOL,1);
fprintf('%i outliers removed (%1.1f%%), %i remain\n',origSz-fnSize,100*((origSz-fnSize)/origSz),fnSize);

%% plot all seasons with mast cond
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
close all;
mastCond = {[0,1],0,1};
mastTitle = {'All Years','Non-mast','Mast'};
doSave = 0;
for iMast = 1:3
    ff(500,900);
    rows = 5;
    cols = 2;
    for iSeason = 1:5
        if iSeason == 1
            useIds = find(T_SOL.isSunrise==1 & ismember(T_SOL.is_mast,mastCond{iMast}));
            SOLs_sunrise = T_SOL.SOL(useIds);
            sunrise_means = mean([T_SOL.awakeIdx(useIds),T_SOL.asleepIdx(useIds)],2);

            useIds = find(T_SOL.isSunrise==0 & ismember(T_SOL.is_mast,mastCond{iMast}));
            SOLs_sunset = T_SOL.SOL(useIds);
            sunset_means = mean([T_SOL.awakeIdx(useIds),T_SOL.asleepIdx(useIds)],2);
            seasonLabel = 'All Seasons';
            sunriseColor = repmat(0.8,[1,3]);
        else
            useIds = find(T_SOL.isSunrise==1 & T_SOL.season == iSeason - 1 & ismember(T_SOL.is_mast,mastCond{iMast}));
            SOLs_sunrise = T_SOL.SOL(useIds);
            sunrise_means = mean([T_SOL.awakeIdx(useIds),T_SOL.asleepIdx(useIds)],2);

            useIds = find(T_SOL.isSunrise==0 & T_SOL.season == iSeason - 1 & ismember(T_SOL.is_mast,mastCond{iMast}));
            SOLs_sunset = T_SOL.SOL(useIds);
            sunset_means = mean([T_SOL.awakeIdx(useIds),T_SOL.asleepIdx(useIds)],2);
            seasonLabel = seasonLabels{iSeason-1};
            sunriseColor = colors(iSeason-1,:);
        end

        subplot(rows,cols,prc(cols,[iSeason,1]));
        binEdges = linspace(0,250,50);
        histogram(SOLs_sunrise,binEdges,'FaceColor',sunriseColor);
        hold on;
        histogram(SOLs_sunset,binEdges,'FaceColor','k');
        set(gca,'fontsize',14);
        if iSeason == 5
            xlabel('Duration (min.)');
        end
        ylabel('Frequency');
        if iSeason == 1
            title({'Latency',sprintf('%s-%s',seasonLabel,mastTitle{iMast})});
        else
            title(sprintf('%s-%s',seasonLabel,mastTitle{iMast}));
        end
        grid on;
        legend({'Sunrise','Sunset'},'fontsize',11);
        legend boxoff;

        subplot(rows,cols,prc(cols,[iSeason,2]));
        binEdges = linspace(720-400,720+400,50);
        histogram(sunrise_means,binEdges,'FaceColor',sunriseColor);
        hold on;
        histogram(sunset_means,binEdges,'FaceColor','k');
        set(gca,'fontsize',14);
        if iSeason == 5
            xlabel('Rel. Time (min.)');
        end
    % %     ylabel('Frequency');
        xticks([720-300,720,720+300]);
        xticklabels({'-300','0','+300'});
        if iSeason == 1
            title({'Onset/Offset',sprintf('%s-%s',seasonLabel,mastTitle{iMast})});
        else
            title(sprintf('%s-%s',seasonLabel,mastTitle{iMast}));
        end
        grid on;
    % %     legend({'QB-AB','AB-QB'},'fontsize',11,'autoupdate','off');
    % %     legend boxoff;
        xline(720,'k-');

        drawnow;
    end

    if doSave
        print(gcf,'-painters','-depsc',fullfile(exportPath,sprintf('%s-%s.eps','T_SOL_cleanHistograms',mastTitle{iMast}))); % required for vector lines
        saveas(gcf,fullfile(exportPath,sprintf('%s-%s.jpg','T_SOL_cleanHistograms',mastTitle{iMast})),'jpg');
        close(gcf);
    end
end
%% T_SOL table and p-value matrix (3/3)
rowNames = {'Latency to AB','AB Rel. to Sunrise','Latency to QB','QB Rel. to Sunset'};
varTypes = {'string','string','string','string','string','string'};
T_SOL_summary = table('Size',[4,6],'VariableNames',{'Description',seasonLabels{:},'All'},'VariableType',varTypes);
isSunriseArr = [1,0];
iRow = 0;
for ii = 1:numel(rowNames)
    T_SOL_summary.Description(ii) = rowNames{ii}; 
end
for iSeason = 1:5
    for iSun = 1:2
        if iSeason == 5
            useRows = find(T_SOL.isSunrise == isSunriseArr(iSun));
        else
            useRows = find(T_SOL.season == iSeason & T_SOL.isSunrise == isSunriseArr(iSun));
        end
        asleepAwakeMean = mean([T_SOL.awakeIdx(useRows),T_SOL.asleepIdx(useRows)],2);
        SOLs = T_SOL.SOL(useRows);
        
        if iSun == 1
            T_SOL_summary(1,iSeason+1) = {sprintf('%1.2f ± %1.2f',mean(SOLs),std(SOLs))};
            T_SOL_summary(2,iSeason+1) = {sprintf('%1.2f ± %1.2f',mean(asleepAwakeMean)-720,std(asleepAwakeMean))};
        else
            T_SOL_summary(3,iSeason+1) = {sprintf('%1.2f ± %1.2f',mean(SOLs),std(SOLs))};
            T_SOL_summary(4,iSeason+1) = {sprintf('%1.2f ± %1.2f',mean(asleepAwakeMean)-720,std(asleepAwakeMean))};
        end
    end
end
writetable(T_SOL_summary,'T_SOL_summary.xlsx');

cmap = parula(1000);
cmap = [cmap(100:900,:);0,0,0];
clc
seasonAbbr = {'Wi','Sp','Su','Au'};
close all;
ff(600,120);
for iSubplot = 1:4
    pMat_sol = NaN(4,4);
    for iSeason = 1:4
        for kSeason = iSeason:4
            if iSubplot == 1
                s1Ids = find(T_SOL.season == iSeason & T_SOL.isSunrise == 1);
                s2Ids = find(T_SOL.season == kSeason & T_SOL.isSunrise == 1);
                y = [T_SOL.SOL(s1Ids);T_SOL.SOL(s2Ids)];
            elseif iSubplot == 2
                s1Ids = find(T_SOL.season == iSeason & T_SOL.isSunrise == 1);
                s2Ids = find(T_SOL.season == kSeason & T_SOL.isSunrise == 1);
                y = [mean([T_SOL.awakeIdx(s1Ids),T_SOL.asleepIdx(s1Ids)],2);mean([T_SOL.awakeIdx(s2Ids),T_SOL.asleepIdx(s2Ids)],2)];
            elseif iSubplot == 3
                s1Ids = find(T_SOL.season == iSeason & T_SOL.isSunrise == 0);
                s2Ids = find(T_SOL.season == kSeason & T_SOL.isSunrise == 0);
                y = [T_SOL.SOL(s1Ids);T_SOL.SOL(s2Ids)];
            else
                s1Ids = find(T_SOL.season == iSeason & T_SOL.isSunrise == 0);
                s2Ids = find(T_SOL.season == kSeason & T_SOL.isSunrise == 0);
                y = [mean([T_SOL.awakeIdx(s1Ids),T_SOL.asleepIdx(s1Ids)],2);mean([T_SOL.awakeIdx(s2Ids),T_SOL.asleepIdx(s2Ids)],2)];
            end
            group = [zeros(size(s1Ids));ones(size(s2Ids))];
            pMat_sol(iSeason,kSeason) = anova1(y,group,'off');
        end
    end
    disp(rowNames{iSubplot});
    flip(pMat_sol)
    writematrix(flip(pMat_sol),fullfile(exportPath,sprintf('%s.csv',rowNames{iSubplot})));
    subplot(1,4,iSubplot);
    alphaData = ~isnan(pMat_sol);
    imagesc(pMat_sol,'AlphaData',alphaData);
    xticks(1:4);
    yticks(xticks);
    xticklabels(seasonAbbr);
    yticklabels(seasonAbbr);
    title(rowNames{iSubplot});
    colormap(cmap); %flip(magma)
    caxis([0 0.0499]);
    set(gca,'fontsize',12);
    set(gca,'ydir','normal');
    drawnow;
end
cb = cbAside(gca,'p-value','k');
cb.FontSize = 12;
cb.TickLabels = [0 0.05];

doSave = 1;
if doSave
    print(gcf,'-painters','-depsc',fullfile(exportPath,'T_SOL_summary_pMatrix.eps')); % required for vector lines
    saveas(gcf,fullfile(exportPath,'T_SOL_summary_pMatrix.jpg'),'jpg');
    close(gcf);
end

%% (1/3) Probability Density histograms
doSave = 1;
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
nBins_sd = 9; % mins
showMinutes = 60;
sleepThresh = 5; % binWidth below
ssBeforeAfterPad = 30; % minutes
maxY = 0.25;
% main plot
close all
h = ff(450,250);
lns = [];
seasonLabels = {'Winter','Spring','Summer','Autumn'};
sleepDurations_season = {};
countsMat = [];

for iSeason = 1:4
    uniqueSqs = unique(trans_is);
    sleepDurations = [];
    transStarts = [];
    for iSq = uniqueSqs
        theseSleepTrans = find(trans_is==iSq & trans_to==0 & ismember(trans_on,useDoys{iSeason}));
        
        % mast: [2014,2019], nmast: [2015:2018,2020]
        %         theseSleepTrans = find(trans_is==iSq & trans_to==0 & ismember(trans_on,useDoys{iSeason})...
        %             & ismember(trans_yr,[2014,2019]));
        
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
            patch([sleepThresh,max(binEdges_sd),max(binEdges_sd),sleepThresh],[0,0,1,1],'black','FaceAlpha',.1,'EdgeColor','none');
            hold on;
            xline(sleepThresh,'k:');
            text(sleepThresh,maxY-0.02,strcat(sprintf('%i minutes',sleepThresh),'\rightarrow'),'horizontalalignment','right','fontsize',14);
        end
        histogram(sleepDurations,binEdges_sd,'Normalization','probability','EdgeColor',colors(iSeason,:),'DisplayStyle','Stairs','lineWidth',4,'EdgeAlpha',0.75);
        lns(iSeason) = plot(-1,-1,'-','lineWidth',3,'color',colors(iSeason,:));
        counts = histcounts(sleepDurations,binEdges_sd,'Normalization','probability');
        countsMat(iSeason,:) = counts;
        
        % find center of 'sleep' peak
% % % %         overSample = 10000;
% % % %         countsSmooth = smoothdata(equalVectors(counts,overSample),'gaussian',overSample/10);
% % % %         binsSmooth = equalVectors(binEdges,overSample);
% % % %         locs = peakseek(countsSmooth);
% % % %         maxBin = binsSmooth(locs(end));
% % % %         text(17,maxY-iSeason*0.015,strcat(sprintf('%s: %1.2f',seasonLabels{iSeason}(1:2),maxBin),'\rightarrow'),'horizontalalignment','right','fontsize',12,'color',colors(iSeason,:));
        
        %         xline(median(sleepDurations),':','lineWidth',4,'color',colors(iSeason,:));
% % % %         xline(maxBin,':','lineWidth',2,'color',colors(iSeason,:));
        set(gca,'xscale','log');
        drawnow;
    end
end

set(gca,'fontsize',14);
xlabel('QB Duration (minutes)');
ylabel('Probability');
ylim([0 maxY]);
title('Dark QB Duration');
grid on;
legend(lns,seasonLabels);
legend box off

if doSave
    print(gcf,'-painters','-depsc',fullfile(exportPath,'DarkQBDuration.eps')); % required for vector lines
    saveas(gcf,fullfile(exportPath,'DarkQBDuration.jpg'),'jpg');
    close(h);
end

%% (2/3) setup pMat
nSurr = 1000;
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
%% (3/3) plot pMat
doSave = 1;
cmap = parula(1000);
cmap = [cmap(100:900,:);0,0,0];
fs = 12;
close all;
h = ff(1000,100);
seasonAbbr = {'Wi','Sp','Su','Au'};
for iBin = 1:size(pMat_sd,3)
    subplot(1,size(pMat_sd,3),iBin);
    ps = squeeze(pMat_sd(:,:,iBin));
    imagesc(ps,'AlphaData',~isnan(ps));
    caxis([0,0.0499]);
    colormap(cmap); % flip(magma)
    title(sprintf('%1.0f-%1.0f min.',binEdges_sd(iBin),binEdges_sd(iBin+1)));
    xticks(1:4);
    yticks(1:4);
    xticklabels(seasonAbbr);
    yticklabels(seasonAbbr);
    set(gca,'fontsize',fs);
    set(gca,'ydir','normal')
end
cb = cbAside(gca,'p-value','k');
cb.TickLabels = [0 0.05];
cb.FontSize = fs;
if doSave
    print(gcf,'-painters','-depsc',fullfile(exportPath,'DarkQBDuration_pMat.eps')); % required for vector lines
    saveas(gcf,fullfile(exportPath,'DarkQBDuration_pMat.jpg'),'jpg');
    close(h);
end

%% MEAN AND PROBABILITY QB HEATMAP
doSave = 1;
windowSize = 7;
binWidth = 5; % minutes
bins = round(linspace(1,366,60));
yearDoys = 1:366;
nBins = 144;

seasonMeanDoys = [];
for iSeason = 1:4
    seasonMeanDoys(iSeason) = useDoys{iSeason}(round(numel(useDoys{iSeason})/2));
end
allSunrise = [];
allSunset = [];
for iBin = 1:numel(bins)
    allSunrise(iBin) = 0 - (secDay(Tss.noon(bins(iBin))) - secDay(Tss.sunrise(bins(iBin)))) ./ (3600); % if nBins=144
    allSunset(iBin) = (secDay(Tss.sunset(bins(iBin))) - secDay(Tss.noon(bins(iBin)))) ./ (3600);
end
seasonTicks = {'Winter','Spring','Summer','Autumn'};

close all
h = ff(800,600);
for iPlot = 1:2
    allTransHist = NaN(numel(bins),nBins);
    for iBin = 1:numel(bins)
        theseDoys = circshift(1:366,windowSize-bins(iBin));
        theseDoys = theseDoys(1:windowSize*2);
        transStarts = [];
        iDoyCount = 0;
        histDoyArr = [];
        
        theseSleepTrans = find(trans_to==0 & ismember(trans_on,theseDoys)); %  & ~ismember(trans_yr,[2014,2019])
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
            if iPlot == 1
                transStarts = theseTransStart(theseDurations > 0);
            else
                transStarts = theseTransStart(theseDurations > binWidth);
            end
            [counts,binEdges,C] = histcounts(transStarts,linspace(0,86400,nBins+3),'normalization','probability'); %'BinMethod','sqrt'
            
            % use bin indexes to get mean sleep duration, overwrite counts
            if iPlot == 1
                %                 theseDurations = theseDurations(theseDurations > binWidth);
                for ii = 1:numel(counts)
                    counts(ii) = median(rmoutliers(theseDurations(C == ii)));
                end
            end
            
            counts = counts(2:end-1); % rm edges
            binShift = round(((86400/2) - mean(secDay(Tss.noon(theseDoys)))) / (86400/nBins));
            allTransHist(iBin,:) = circshift(counts,binShift);
        end
        disp(iBin);
    end
    
    subplot(2,1,iPlot);
    xSolarNoon = linspace(-12,12,size(allTransHist,2));
    pcolor(xSolarNoon,bins,allTransHist);
    shading interp;
    xticks(-12:2:12);
    useyticks = round(linspace(1,366,16));
    yticks(useyticks);
    
    colormap(magma(6));
    set(gca,'ydir','normal');
    
    useyticklabels = compose("%i",useyticks);
    for iSeason = 1:4
        targetDoy = useDoys{iSeason}(round(numel(useDoys{iSeason})/2));
        yidx = closest(yticks,targetDoy);
        useyticklabels{yidx} = convertStringsToChars("\textbf{" + seasonTicks{iSeason} + "}");
    end
    yticklabels(useyticklabels);
    set(gca,'fontsize',14,'TickLabelInterpreter','latex','FontName','phv');
    xlabel("Hours Relative to Solar Noon");
    ylabel('Julian Day of Year');
    c = colorbar;
    
    if iPlot == 1
        title('Median QB Duration');
        ylabel(c,'Minutes','fontsize',14);
        %         caxis([0 25]);
        caxisauto(allTransHist,1)
    else
        title(sprintf('Probability of QB >%i min.',binWidth));
        ylabel(c,'Probability','fontsize',14);
        %         caxis([0.002 .014]);
        caxisauto(allTransHist,1)
    end
    
    %     caxis(caxis*0.8) % amplify colors
    text(0,290,'No Data','verticalalignment','middle','horizontalalignment','center','fontsize',14);
    
    hold on;
    plot(allSunrise,bins,'color',[1 1 1 0.6],'linewidth',4);
    plot(allSunset,bins,'color',[1 1 1 0.6],'linewidth',4);
    hold off;
end

if doSave
    print(gcf,'-painters','-depsc',fullfile(exportPath,'QBTransitions.eps')); % required for vector lines
    saveas(gcf,fullfile(exportPath,'QBTransitions.jpg'),'jpg');
    close(gcf);
end

%% this is all the transitions (x-y plot) for all seasons
% sunrise/set centered on top, solar noon on bottom
Tss = makeTss(2014:2020);
months =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

close all
ff(750,600);
binEdges = linspace(1,86400,100); % rm 0 entries, are those subsequent animal entries?

nHalfWindow = 30;
allDoys = 1:366;
colors = seasonColors(1:366); % this will use all doys because it's windowed
op = 1;
nS = 1;
t = linspace(0,24,numel(binEdges)-1);

useylim = 0.035;
for iSun = 1:2
    for iDoy = 1:366
        shiftDoys = circshift(allDoys,-iDoy+1+nHalfWindow);
        theseDoys = shiftDoys(1:nHalfWindow*2+1);
        shiftBy = closest(t,mean(secDay(Tss.noon(Tss.doy == iDoy)),1)/60/60) + round(numel(t)/2);
        if iSun == 1
            shiftBy = closest(t,mean(secDay(Tss.sunrise(Tss.doy == iDoy)),1)/60/60) + round(numel(t)/2);
        end
        
        subplot(2,2,prc(2,[iSun,1]));
        useIds = trans_to==1 & ismember(trans_on,theseDoys);
        counts = histcounts(trans_at(useIds),binEdges) / sum(useIds);
        plot(t,circshift(imgaussfilt(counts,nS,'padding','circular'),-shiftBy),'color',[colors(iDoy,:),op]);
        hold on;
        if iDoy == 1
            title('Transition to Awake');
            xlim([0 24]);
            ylim([0 useylim]);
            ylabel('Probability')
            if iSun == 1
                xlabel('Relative to Sunrise (hrs)');
            else
                xlabel('Relative to Solar Noon (hrs)');
            end
            xticks(0:6:24);
            xticklabels({'±12','-6','0','+6','±12'});
            set(gca,'fontsize',16)
            grid on;
            c = colorbar('location','southoutside');
            colormap(colors);
            c.Limits = [0,1];
            c.Ticks = linspace(0,1,12);
            c.TickLabels = months;
            c.TickDirection = 'out';
            c.FontSize = 11;
            plot([12 12],ylim,'k:');
            plot([36 36],ylim,'k:');
        end
        
        if iSun == 1
            shiftBy = closest(t,mean(secDay(Tss.sunset(Tss.doy == iDoy)),1)/60/60) + round(numel(t)/2);
        end
        
        subplot(2,2,prc(2,[iSun,2]));
        useIds = trans_to==0 & ismember(trans_on,theseDoys);
        counts = histcounts(trans_at(useIds),binEdges) / sum(useIds);
        plot(t,circshift(imgaussfilt(counts,nS,'padding','circular'),-shiftBy),'color',[colors(iDoy,:),op]);
        hold on;
        if iDoy == 1
            title('Transition to Asleep');
            xlim([0 24]);
            ylim([0 useylim]);
            xticks(0:6:24);
            xticklabels({'±12','-6','0','+6','±12'});
            ylabel('Probability');
            if iSun == 1
                xlabel('Relative to Sunset (hrs)');
            else
                xlabel('Relative to Solar Noon (hrs)');
            end
            set(gca,'fontsize',16)
            grid on;
            c = colorbar('location','southoutside');
            colormap(colors);
            c.Limits = [0,1];
            c.Ticks = linspace(0,1,12);
            c.TickLabels = months;
            c.TickDirection = 'out';
            c.FontSize = 11;
            plot([12 12],ylim,'k:');
            plot([36 36],ylim,'k:');
        end
    end
    drawnow;
end

% alpha has to be set to 1 then changed in the eps, select > same > stroke
% weight, opacity = 15%
doSave = 0;
if doSave
    print(gcf,'-painters','-depsc',fullfile(exportPath,'SleepTransitions.eps')); % required for vector lines
    saveas(gcf,fullfile(exportPath,'SleepTransitions.jpg'),'jpg');
    close(gcf);
end

%% Probability Density plots OLD
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);

showHours = 3; % hours
binWidth = 10; % mins
% main plot
close all
h = ff(1400,350);
lw1 = 3;
lw2 = 1;
rows = 2;
cols = 4;
lns = [];
seasonLabels = {'Winter','Spring','Summer','Autumn'};
sleepDurations_season = {};
seasonMarkers = ['*','^','o','+'];
clc
ys = {};
allTransStarts = {};
subplot(131);
for iSeason = 1:4
    uniqueSqs = unique(trans_is);
    sleepDurations = [];
    transStarts = [];
    for iSq = uniqueSqs
        theseSleepTrans = find(trans_is==iSq & trans_to==0 & ismember(trans_on,useDoys{iSeason}));
        
        % mast: [2014,2019], nmast: [2015:2018,2020]
        %         theseSleepTrans = find(trans_is==iSq & trans_to==0 & ismember(trans_on,useDoys{iSeason})...
        %             & ismember(trans_yr,[2014,2019]));
        
        % START light/dark transitions
        % % % %         theseDoys = trans_on(theseSleepTrans);
        % % % %         meanSunrise = mean(secDay(Tss.sunrise(theseDoys)));
        % % % %         meanSunset = mean(secDay(Tss.sunset(theseDoys)));
        % % % %         temp = [];
        % % % %         for ii = theseSleepTrans
        % % % %             if trans_at(ii) > meanSunset || trans_at(ii) < meanSunrise
        % % % %                 temp = [temp ii];
        % % % %             end
        % % % %         end
        % % % %         if isempty(temp)
        % % % %             continue;
        % % % %         end
        % % % %         theseSleepTrans = temp;
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
            % !! only save long transitions for next figure
            transStarts = [transStarts theseTransStart(theseDurations > binWidth)];
            
            sleepDurations = [sleepDurations theseDurations];
        end
    end
    
    if ~isempty(sleepDurations)
        useDurations = sleepDurations(sleepDurations < showHours * 60);
        sleepDurations_season{iSeason} = useDurations;
        allTransStarts{iSeason} = transStarts;
        %         binEdges = linspace(0,showHours*60,nBins);
        binEdges = 0:binWidth:showHours*60;
        [y,~,bin] = histcounts(useDurations,binEdges,'Normalization','pdf');
        x = binEdges(1:end-1) + median(diff(binEdges))/2;
        
        % old markers: seasonMarkers(iSeason)
        lns(iSeason) = plot(x,y,seasonMarkers(iSeason),'color',colors(iSeason,:),'markerSize',10,'lineWidth',3);
        hold on;
        set(gca,'yscale','log');
        
        for iBin = 1:numel(x)
            ys{iSeason,iBin} = useDurations(bin==iBin);
        end
        % % % %         if iSeason == 4
        % % % %             pvals = [];
        % % % %             clc
        % % % %             for iBin = 1:numel(x)
        % % % %                 anova_y = [];
        % % % %                 anova_group = [];
        % % % %                 for anSeason = [2,4]
        % % % %                     anova_y = [anova_y ys{anSeason,iBin}];
        % % % %                     anova_group = [anova_group anSeason*ones(size(ys{anSeason,iBin}))];
        % % % %                 end
        % % % %                 p = anova1(anova_y',anova_group','off');
        % % % %                 fprintf("iBin = %i, p = %1.4f\n",iBin,p);
        % % % %                 pvals(iBin) = p;
        % % % %             end
        % % % %             yyaxis right;
        % % % %             plot(x,pvals,'-','color',repmat(0.7,[3,1]),'lineWidth',3);
        % % % %             hold on;
        % % % %
        % % % %             pvals = [];
        % % % %             clc
        % % % %             for iBin = 1:numel(x)
        % % % %                 anova_y = [];
        % % % %                 anova_group = [];
        % % % %                 for anSeason = [1,3]
        % % % %                     anova_y = [anova_y ys{anSeason,iBin}];
        % % % %                     anova_group = [anova_group anSeason*ones(size(ys{anSeason,iBin}))];
        % % % %                 end
        % % % %                 p = anova1(anova_y',anova_group','off');
        % % % %                 fprintf("iBin = %i, p = %1.4f\n",iBin,p);
        % % % %                 pvals(iBin) = p;
        % % % %             end
        % % % %             yyaxis right;
        % % % %             plot(x,pvals,'k-','lineWidth',3);
        % % % % %             plot(x,pvals,'--','color',[colors(3,:),0.5],'lineWidth',3);
        % % % %             set(gca,'ycolor','k');
        % % % %             ylim([0 1]);
        % % % %             ylabel('p-value');
        % % % %             yyaxis left;
        % % % %         end
        
        onlyUseIds = find(x >= 30);
        coeffs = polyfit(x(onlyUseIds),log(y(onlyUseIds)),1);
        %         fprintf("%s \tau = %1.2f\n",seasonLabels{iSeason},abs(1/coeffs(1)));
        fittedX = linspace(x(min(onlyUseIds)), x(max(onlyUseIds)), 100);
        fittedLogY = polyval(coeffs, fittedX);
        b = exp(coeffs(2));
        m = exp(coeffs(1));
        % Get y fitted another way
        fittedY2 = b * m .^ fittedX; % b * m .^ fittedX;
        plot(fittedX,fittedY2,'-','color',colors(iSeason,:),'lineWidth',1);
    end
end

% do stats
% % % % ylab = flip(linspace(0.4,0.8,6));
% % % % statCount = 0;
% % % % fs = 12;
% % % % statColor = repmat(0.5,[1,4]);
% % % % for ii = 1:4
% % % %     for jj = ii+1:4
% % % %         % if p < 0.05 dist. are different
% % % %         [h,p] = kstest2(sleepDurations_season{ii},sleepDurations_season{jj});
% % % %         sigStr = '';
% % % %         if p < 0.05 && p >= 0.01
% % % %             sigStr = '*';
% % % %         elseif p < 0.01 && p >= 0.001
% % % %             sigStr = '**';
% % % %         elseif p < 0.001
% % % %             sigStr = '***';
% % % %         end
% % % %         statCount = statCount + 1;
% % % %         sprintfString = sprintf('%s-%s %sp = %1.2e',seasonLabels{ii}(1:3),seasonLabels{jj}(1:3),sigStr,p);
% % % %         disp(sprintfString);
% % % % %         text(4,ylab(statCount),sprintfString,'fontsize',fs,'color',statColor);
% % % %     end
% % % % end

legend(lns,seasonLabels);
xlim([0 showHours*60 + median(diff(binEdges))]);
xticks(x);
set(gca,'fontsize',14);
xlabel('Duration t (min)');
ylabel('Probability Density P');
grid off; grid on;
title('QB Periods');
ylimVals = ylim;
xlimVals = xlim;

% % doSave = 1;
% % if doSave
% %     print(gcf,'-painters','-depsc',fullfile(exportPath,'ProbabilityDensity.eps')); % required for vector lines
% %     saveas(gcf,fullfile(exportPath,'ProbabilityDensity.jpg'),'jpg');
% %     close(gcf);
% % end

%% OLD LINE PLOTS
% transition times +-12 from solar noon
subplot(122);
shiftBy = [];
% circshift(imgaussfilt(counts,nS,'padding','circular'),-shiftBy)
% shiftBy = closest(t,mean(secDay(Tss.sunrise(Tss.doy == iDoy)),1)/60/60) + round(numel(t)/2);
% close all
% ff(400,300);
for iSeason = 1:4
    if isempty(allTransStarts{iSeason})
        continue;
    end
    shiftBy(iSeason) = round(mean(secDay(Tss.noon(useDoys{iSeason}))));
    shiftedTrans = allTransStarts{iSeason}-shiftBy(iSeason);
    %     shiftedBins = linspace(min(shiftedTrans),max(shiftedTrans),nBins);
    [counts,binEdges] = histcounts(shiftedTrans,'normalization','probability','BinMethod','sqrt');
    counts = imgaussfilt(counts(2:end-1),1,'padding','circular');
    %     circCounts = circshift(counts,round((-shiftBy(iSeason)/1440)*100)); % imgaussfilt(counts,1,'padding','circular'),
    %     circx = linspace(-pi,pi,numel(counts));
    %     polarplot(circx,counts,'linewidth',2,'color',colors(iSeason,:));
    x = linspace(-12,12,numel(counts));
    plot(x,counts,'-','color',colors(iSeason,:),'linewidth',3);
    hold on;
end
set(gca,'fontsize',14);
title("AB\rightarrowQB Transitions");
xlabel("Hours to Solar Noon");
xticks(-12:6:12);
xlim([-12 12]);
ylabel('Probability');
grid on;

if doSave
    print(gcf,'-painters','-depsc',fullfile(exportPath,'QBTransitions.eps')); % required for vector lines
    saveas(gcf,fullfile(exportPath,'QBTransitions.jpg'),'jpg');
    close(gcf);
end

%% mast version
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);

% showHours = 3;
% binWidth = 20;
h = ff(500,500);
all_x = {};
all_y = {};
lnmast = [];
for iSeason = 1:4
    subplot(2,2,iSeason);
    lns = NaN(2,1);
    sleepDurations_mast = {};
    for iMast = 1:2
        if iMast == 1
            useYears = [2015:2018,2020];
            op = 1;
        else
            useYears = [2014,2019];
            op = 0.4;
        end
        uniqueSqs = unique(trans_is);
        sleepDurations = [];
        for iSq = uniqueSqs
            theseSleepTrans = find(trans_is==iSq & trans_to==0 &...
                ismember(trans_on,useDoys{iSeason}) & ismember(trans_yr,useYears));
            
            % test light/dark transitions
            theseDoys = trans_on(theseSleepTrans);
            meanSunrise = mean(secDay(Tss.sunrise(theseDoys)));
            meanSunset = mean(secDay(Tss.sunset(theseDoys)));
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
            
            allEntries = find(trans_is==iSq);
            % must end with awake transition to get time
            if trans_to(allEntries(end)) == 0
                theseSleepTrans = theseSleepTrans(1:end-1);
            end
            
            theseTransStart = trans_at(theseSleepTrans); % asleep
            theseTransEnd = trans_at(theseSleepTrans+1); % awake
            theseDurations = (theseTransEnd - theseTransStart) / 60;
            % fix transitions around next day (in seconds)
            for ii = find(theseDurations < 0)
                theseDurations(ii) = theseDurations(ii) + 1440;
            end
            sleepDurations = [sleepDurations theseDurations];
        end
        useDurations = sleepDurations(sleepDurations < showHours * 60);
        
        binEdges = 0:binWidth:showHours*60;
        [y,~,bin] = histcounts(useDurations,binEdges,'Normalization','pdf');
        x = binEdges(1:end-1) + median(diff(binEdges))/2;
        
        sleepDurations_mast{iMast} = useDurations;
        sleepBins_mast{iMast} = bin;
        
        all_x{iSeason,iMast} = x;
        all_y{iSeason,iMast} = y;
        
        if iMast == 1
            lns(iSeason) = plot(x,y,'.','color',colors(iSeason,:),'markerSize',30);
        else
            lns(iSeason) = plot(x,y,'x','color',colors(iSeason,:),'markerSize',12,'linewidth',3);
            if iSeason == 1
                continue;
            end
            for iBin = 1:numel(x)
                nmastDurs = sleepDurations_mast{1};
                nmastDurs = nmastDurs(sleepBins_mast{1}==iBin);
                mastDurs = sleepDurations_mast{2};
                mastDurs = mastDurs(sleepBins_mast{2}==iBin);
                group = [zeros(numel(nmastDurs),1);ones(numel(mastDurs),1)];
                p = anova1([nmastDurs,mastDurs]',group,'off');
                fprintf("iSeason: %i, iBin: %i, p = %1.4f\n",iSeason,iBin,p);
                if p < 0.05
                    plot(x(iBin),max(ylim),'r*');
                end
                disp("");
            end
        end
        hold on;
        lnmast(1) = plot(0,0,'k-','lineWidth',4);
        lnmast(2) = plot(0,0,'k:','lineWidth',4);
        hold on;
        grid on;
        set(gca,'yscale','log');
        xlim([0,showHours*60]);
        ylim([10^-5 1])
        %         xlim(xlimVals);
        %         ylim(ylimVals);
    end
end
% legend(lnmast,{'Non-mast','Mast'},'location','southwest');
% set(gca,'fontsize',14);
% xticklabels({});
% yticklabels({});
% grid off;
% box off;

doSave = 1;
if doSave
    print(gcf,'-painters','-depsc',fullfile(exportPath,'ProbabilityDensity_mast.eps')); % required for vector lines
    saveas(gcf,fullfile(exportPath,'ProbabilityDensity_mast.jpg'),'jpg');
    close(gcf);
end

%%
close all
ff(800,200);
for iSeason = 2:4
    subplot(1,3,iSeason-1);
    vals = all_y{iSeason,2} - all_y{iSeason,1};
    binVals = (vals > 0) + ((vals <= 0) * -1);
    y = [sum(vals > 0),sum(vals <= 0)];
    plot(all_x{iSeason,2},binVals,'.','color',colors(iSeason,:),'markersize',30);
    hold on;
    plot(all_x{iSeason,2},binVals,'-','color',[colors(iSeason,:),0.15]);
    ylim([-3 3]);
    yticks([-1,1]);
    yticklabels({'Non-mast','Mast'});
    title(sprintf("%1.2f%% Mast",100*y(1)/sum(y)));
    set(gca,'fontsize',14);
    xlabel('Duration t (min)');
    %     pie(y);
    hold on;
    %     legend({'Mast','Non-mast'});
end
saveas(gcf,fullfile(exportPath,'BinaryQuiescence_mast.jpg'),'jpg');

%% older CDF plots
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);

% main plot
close all
ff(1000,400);
lw1 = 3;
lw2 = 1;
rows = 2;
cols = 4;
subplot(rows,cols,[1 2 5 6]);
lns = [];
showHours = 8;
seasonLabels = {'Winter','Spring','Summer','Autumn'};
sleepDurations_season = {};
for iSeason = 1:4
    uniqueSqs = unique(trans_is);
    sleepDurations = [];
    for iSq = uniqueSqs
        theseSleepTrans = find(trans_is==iSq & trans_to==0 & ismember(trans_on,useDoys{iSeason}));
        if ~isempty(theseSleepTrans)
            
            % test light/dark transitions
            % % % %             theseDoys = trans_on(theseSleepTrans);
            % % % %             meanSunrise = mean(secDay(Tss.sunrise(theseDoys)));
            % % % %             meanSunset = mean(secDay(Tss.sunset(theseDoys)));
            % % % %             temp = [];
            % % % %             for ii = theseSleepTrans
            % % % %                 if trans_at(ii) > meanSunset || trans_at(ii) < meanSunrise
            % % % %                     temp = [temp ii];
            % % % %                 end
            % % % %             end
            % % % %             if isempty(temp)
            % % % %                 continue;
            % % % %             end
            % % % %             theseSleepTrans = temp;
            
            allEntries = find(trans_is==iSq);
            % must end with awake transition to get time
            if trans_to(allEntries(end)) == 0
                theseSleepTrans = theseSleepTrans(1:end-1);
            end
            
            theseTransStart = trans_at(theseSleepTrans); % asleep
            theseTransEnd = trans_at(theseSleepTrans+1); % awake
            % make minutes, that's the highest resolution from
            % predict_awake anyways
            theseDurations = (theseTransEnd - theseTransStart) / 60;
            % fix transitions around next day (in seconds)
            for ii = find(theseDurations < 0)
                theseDurations(ii) = theseDurations(ii) + 1440;
            end
            sleepDurations = [sleepDurations theseDurations];
        end
    end
    
    % uncomment for PDF
    %     [counts,edges] = histcounts(sleepDurations/60,linspace(0,8,(8/24)*1440));
    %     plot(edges(1:end-1),normalize(smoothdata(counts,'gaussian',5),'range'),'color',colors(iSeason,:),'linewidth',2);
    %     hold on;
    
    % uncomment for CDF
    sleepDurations_season{iSeason} = sleepDurations;
    if ~isempty(sleepDurations)
        [Fout,x,Flo,Fup] = ecdf(sleepDurations/60); % hours
        lns(iSeason) = plot(x,Fout,'color',colors(iSeason,:),'linewidth',lw1);
        hold on;
        plot(x,Flo,':','color',colors(iSeason,:),'linewidth',lw2);
        plot(x,Fup,':','color',colors(iSeason,:),'linewidth',lw2);
    end
end

% do stats
ylab = flip(linspace(0.4,0.8,6));
statCount = 0;
fs = 12;
statColor = repmat(0.5,[1,4]);
for ii = 1:4
    for jj = ii+1:4
        % if p < 0.05 dist. are different
        [h,p] = kstest2(sleepDurations_season{ii},sleepDurations_season{jj});
        sigStr = '';
        if p < 0.05 && p >= 0.01
            sigStr = '*';
        elseif p < 0.01 && p >= 0.001
            sigStr = '**';
        elseif p < 0.001
            sigStr = '***';
        end
        statCount = statCount + 1;
        text(4,ylab(statCount),sprintf('%s-%s %sp = %1.2e',seasonLabels{ii}(1:3),...
            seasonLabels{jj}(1:3),sigStr,p),'fontsize',fs,'color',statColor);
    end
end


xlim([0 showHours]);
set(gca,'fontsize',16);
xlabel('Asleep Length (hrs)');
ylabel('Probability');
grid on;
title('Asleep Periods');
legend(lns,seasonLabels,'location','southeast','AutoUpdate','off');

% mast plots
useSubplots = [3,4,7,8];
all_sleepDurations = [];
for iSeason = 1:4
    lns = NaN(2,1);
    sleepDurations_mast = {};
    for iMast = 1:2
        if iMast == 1
            useYears = [2015:2020];
            op = 1;
        else
            useYears = [2014,2019];
            op = 0.4;
        end
        subplot(rows,cols,useSubplots(iSeason));
        uniqueSqs = unique(trans_is);
        sleepDurations = [];
        for iSq = uniqueSqs
            theseSleepTrans = find(trans_is==iSq & trans_to==0 &...
                ismember(trans_on,useDoys{iSeason}) & ismember(trans_yr,useYears));
            
            % test light/dark transitions
            % % % %             theseDoys = trans_on(theseSleepTrans);
            % % % %             meanSunrise = mean(secDay(Tss.sunrise(theseDoys)));
            % % % %             meanSunset = mean(secDay(Tss.sunset(theseDoys)));
            % % % %             temp = [];
            % % % %             for ii = theseSleepTrans
            % % % %                 if trans_at(ii) > meanSunset || trans_at(ii) < meanSunrise
            % % % %                     temp = [temp ii];
            % % % %                 end
            % % % %             end
            % % % %             if isempty(temp)
            % % % %                 continue;
            % % % %             end
            % % % %             theseSleepTrans = temp;
            
            allEntries = find(trans_is==iSq);
            % must end with awake transition to get time
            if trans_to(allEntries(end)) == 0
                theseSleepTrans = theseSleepTrans(1:end-1);
            end
            
            theseTransStart = trans_at(theseSleepTrans); % asleep
            theseTransEnd = trans_at(theseSleepTrans+1); % awake
            theseDurations = (theseTransEnd - theseTransStart) / 60;
            % fix transitions around next day (in seconds)
            for ii = find(theseDurations < 0)
                theseDurations(ii) = theseDurations(ii) + 1440;
            end
            sleepDurations = [sleepDurations theseDurations];
        end
        all_sleepDurations = [all_sleepDurations sleepDurations];
        
        
        sleepDurations_mast{iMast} = sleepDurations;
        if ~isempty(sleepDurations)
            [Fout,x,Flo,Fup] = ecdf(sleepDurations/60);
            lns(iMast) = plot(x,Fout,'color',[colors(iSeason,:),op],'linewidth',lw1);
            hold on;
            plot(x,Flo,':','color',[colors(iSeason,:),op],'linewidth',lw2);
            plot(x,Fup,':','color',[colors(iSeason,:),op],'linewidth',lw2);
        end
        xlim([0 showHours]);
        ylim([0 1]);
        set(gca,'fontsize',12);
        xlabel('Asleep Length (hrs)');
        grid on;
        title(seasonLabels{iSeason});
        if ~isnan(lns(1)) && ~isnan(lns(2))
            legend(lns,{'Non-mast','Mast'},'location','southeast','AutoUpdate','off');
        elseif ~isnan(lns(1))
            legend(lns(1),{'Non-mast'},'location','southeast','AutoUpdate','off');
        else
            legend(lns(2),{'Mast'},'location','southeast','AutoUpdate','off');
        end
    end
    
    if ~isempty(sleepDurations_mast{1}) && ~isempty(sleepDurations_mast{2})
        [h,p] = kstest2(sleepDurations_mast{1},sleepDurations_mast{2});
        sigStr = '';
        if p < 0.05 && p >= 0.01
            sigStr = '*';
        elseif p < 0.01 && p >= 0.001
            sigStr = '**';
        elseif p < 0.001
            sigStr = '***';
        end
        
        text(2,mean(ylab),sprintf('%sp = %1.2e',sigStr,p),'fontsize',fs,'color',statColor);
    else
        text(2,mean(ylab),'no mast data','fontsize',fs,'color',statColor);
    end
end

%% test each year
seasonLabels = {'Winter','Spring','Summer','Autumn'};
lw1 = 3;
lw2 = 1;
close all
ff(1200,900);
% iSeason = 4;
for iSeason = 1:4
    subplot(2,2,iSeason);
    lns = [];
    legendLabels = {};
    for iYear = 2014:2020
        if ismember(iYear,[2014,2019])
            useColor = 'r';
        else
            useColor = 'k';
        end
        uniqueSqs = unique(trans_is);
        sleepDurations = [];
        for iSq = uniqueSqs
            theseSleepTrans = find(trans_is==iSq & trans_to==0 &...
                ismember(trans_on,useDoys{iSeason}) & trans_yr == iYear);
            nSq = numel(unique(trans_is(ismember(trans_on,useDoys{iSeason}) & trans_yr == iYear)));
            
            allEntries = find(trans_is==iSq);
            % must end with awake transition to get time
            if trans_to(allEntries(end)) == 0
                theseSleepTrans = theseSleepTrans(1:end-1);
            end
            
            theseTransStart = trans_at(theseSleepTrans); % asleep
            theseTransEnd = trans_at(theseSleepTrans+1); % awake
            theseDurations = (theseTransEnd - theseTransStart) / 60;
            % fix transitions around next day (in seconds)
            for ii = find(theseDurations < 0)
                theseDurations(ii) = theseDurations(ii) + 1440;
            end
            sleepDurations = [sleepDurations theseDurations];
        end
        % uncomment for PDF
        %     [counts,edges] = histcounts(transTimes/3600,linspace(0,8,(8/24)*1440));
        %     plot(edges(1:end-1),normalize(smoothdata(counts,'gaussian',20),'range'),'color',colors(iSeason,:),'linewidth',2);
        %     hold on;
        
        % uncomment for CDF
        if ~isempty(sleepDurations)
            legendLabels{numel(legendLabels)+1} = sprintf('%i, n = %i',iYear,nSq);
            [Fout,x,Flo,Fup] = ecdf(sleepDurations/60);
            lns(numel(lns)+1) = plot(x,Fout,'color',useColor,'linewidth',lw1);
            hold on;
            plot(x,Flo,':','color',useColor,'linewidth',lw2);
            plot(x,Fup,':','color',useColor,'linewidth',lw2);
        end
    end
    legend(lns,legendLabels,'location','southeast');
    xlim([0 8]);
    set(gca,'fontsize',16);
    xlabel('Asleep Length (hrs)');
    ylabel('%');
    grid on;
    title(sprintf('Asleep Periods, %s',seasonLabels{iSeason}));
end