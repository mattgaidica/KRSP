%% !RUNFIG! turn on figure; setup in predict_awake.m
close all
ff(1000,900);
doSave = 1;

cols = 9 + 1 + 8; % =18, col10 is gutter
rows = 3 + 3 + 3 + 1; % =10
figs = {};
figs{1} = [1:9 19:27 37:45];
figs{2} = [55:63 73:81 91:99];
figs{3} = [109:117 127:135 145:153];
figs{4} = 163:171;

figs{5} = [11:14 29:32];
figs{6} = figs{5} + 4;
figs{7} = figs{5} + 36;
figs{8} = figs{6} + 36;
figs{9} = figs{7} + 36;
figs{10} = figs{8} + 36;
figs{11} = figs{9} + 36;
figs{12} = figs{10} + 36;
figs{13} = figs{11} + 36;
figs{14} = figs{12} + 36;

contQB = 5; % minutes
gcaFontSize = 12;
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);

%% !RUNFIG!, set do=1 to recompute, MEAN AND PROBABILITY QB HEATMAP
leftPanelPos = [1 1 1 0.85];
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

% % % % close all
% % % % h = ff(800,600);
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
                theseSleepTrans = find(trans_to==0 & ismember(trans_on,theseDoys) & trans_is==iSq); %  & ~ismember(trans_yr,[2014,2019])
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
            disp(iBin);
        end
    end
end

for iPlot = 1:2
% % % %     subplot(2,1,iPlot);
    subplot(rows,cols,figs{iPlot});
    pos = get(gca,'Position');
    set(gca,'Position',pos.*leftPanelPos + [0 0.05 0 0]);
    
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
    set(gca,'fontsize',gcaFontSize);
    if iPlot == 2
        xlabel("Rel. to Solar Noon (hours)");
    end
    ylabel('Julian Day of Year');
    c = colorbar;
    
    if iPlot == 1
        title('Median QB Duration');
        ylabel(c,'Minutes','fontsize',gcaFontSize);
        %         caxis([0 25]);
        caxisauto(squeeze(allTransHist(iPlot,:,:)),1)
    else
        title(sprintf('Probability of QB >%i min.',contQB));
        ylabel(c,'Probability','fontsize',gcaFontSize);
        %         caxis([0.002 .014]);
        caxisauto(squeeze(allTransHist(iPlot,:,:)),1)
    end
    
    %     caxis(caxis*0.8) % amplify colors
    text(0,300,'No Data','verticalalignment','middle','horizontalalignment','center','fontsize',gcaFontSize);
    
    hold on;
    plot(allSunrise,seasonBins,'color',[1 1 1 0.6],'linewidth',4);
    plot(allSunset,seasonBins,'color',[1 1 1 0.6],'linewidth',4);
    hold off;
end

% % % % if doSave
% % % %     print(gcf,'-painters','-depsc',fullfile(exportPath,'QBTransitions.eps')); % required for vector lines
% % % %     saveas(gcf,fullfile(exportPath,'QBTransitions.jpg'),'jpg');
% % % %     close(gcf);
% % % % end

%% !RUNFIG! (1/3) Probability Density histograms
nBins_sd = 9; % mins
showMinutes = 60;
ssBeforeAfterPad = 30; % minutes
maxY = 0.25;

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
            patch([contQB,max(binEdges_sd),max(binEdges_sd),contQB],[0,0,1,1],'black','FaceAlpha',.1,'EdgeColor','none');
            hold on;
            xline(contQB,'k:');
            text(contQB,maxY-0.02,strcat(sprintf('%i minutes',contQB),'\rightarrow'),'horizontalalignment','right','fontsize',gcaFontSize);
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

set(gca,'fontsize',gcaFontSize);
xlabel('QB Duration (minutes)');
ylabel('Probability');
ylim([0 maxY]);
title('Dark QB Duration');
grid on;
legend(lns,seasonLabels);
legend box off;
hold off;

% % % % if doSave
% % % %     print(gcf,'-painters','-depsc',fullfile(exportPath,'DarkQBDuration.eps')); % required for vector lines
% % % %     saveas(gcf,fullfile(exportPath,'DarkQBDuration.jpg'),'jpg');
% % % %     close(h);
% % % % end

%% (2/3) setup pMat_sd
nSurr = 1000;
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
cmap = parula(1000);
cmap = [cmap(100:900,:);0,0,0];
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
    title(sprintf('%1.0f min',mean([binEdges_sd(iBin),binEdges_sd(iBin+1)])));
    xticks(1:4);
    xtickangle(-90);
    yticks(1:4);
    xticklabels(seasonAbbr);
    yticklabels(seasonAbbr);
    set(gca,'fontsize',gcaFontSize-4);
    set(gca,'ydir','normal')
    pos = get(gca,'Position');
    set(gca,'Position',pos.*[1 1 0.9 0.6]+[0 -0.014 0 0]);
end
cb = cbAside(gca,'p','k');
cb.TickLabels = [0 0.05];
cb.FontSize = gcaFontSize-2;
% % % % if doSave
% % % %     print(gcf,'-painters','-depsc',fullfile(exportPath,'DarkQBDuration_pMat.eps')); % required for vector lines
% % % %     saveas(gcf,fullfile(exportPath,'DarkQBDuration_pMat.jpg'),'jpg');
% % % %     close(h);
% % % % end
%% SOL_T table setup and visualization (1/3), see also: /Users/matt/Documents/MATLAB/KRSP/analyze_circCorrSleep.m
if do
    do = 0;
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
                set(gca,'fontsize',gcaFontSize);
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
                text(max([asleepIdx,awakeIdx])+10,0.5,sprintf('SOL = %i mins',abs(asleepIdx - awakeIdx)),'HorizontalAlignment','left','fontsize',gcaFontSize);
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
end
%% load SOL_T cleaning and seasonal histograms (2/3)
T_SOL = readtable('T_SOL');
origSz = size(T_SOL,1);

T_SOL(T_SOL.awakeIdx == 1,:) = [];
T_SOL(T_SOL.asleepIdx == 1,:) = [];

fnSize = size(T_SOL,1);
fprintf('%i outliers removed (%1.1f%%), %i remain\n',origSz-fnSize,100*((origSz-fnSize)/origSz),fnSize);

%% T_SOL table and p-value matrix (3/3)
doSave = 1;
rowNames = {'Latency to AB','AB Rel. to Sunrise','Latency to QB','QB Rel. to Sunset'};
varTypes = {'string','string','string','string','string','string'};
isSunriseArr = [1,0];
mastCond = {[0,1],0,1};
mastNames = {'All','NonMast','Mast'};
SOLTableFiles = {'T_SOL_summary_All.xlsx','T_SOL_summary_NonMast.xlsx','T_SOL_summary_Mast.xlsx'};
SOLFigureFiles = {'T_SOL_pMatrix_All','T_SOL_pMatrix_NonMast','T_SOL_pMatrix_Mast'};
mastSOL = {};
for iMast = 1:3
    T_SOL_summary = table('Size',[4,6],'VariableNames',[{'Description'},seasonLabels(:)',{'All'}],'VariableType',varTypes);
    for ii = 1:numel(rowNames)
        T_SOL_summary.Description(ii) = rowNames{ii};
    end
    for iSeason = 1:5
        for iSun = 1:2
            if iSeason == 5
                useRows = find(T_SOL.isSunrise == isSunriseArr(iSun) & ismember(T_SOL.is_mast,mastCond{iMast}));
            else
                useRows = find(T_SOL.season == iSeason & T_SOL.isSunrise == isSunriseArr(iSun) & ismember(T_SOL.is_mast,mastCond{iMast}));
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
    writetable(T_SOL_summary,fullfile(exportPath,SOLTableFiles{iMast}));
    
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
                    s1Ids = find(T_SOL.season == iSeason & T_SOL.isSunrise == 1 & ismember(T_SOL.is_mast,mastCond{iMast}));
                    s2Ids = find(T_SOL.season == kSeason & T_SOL.isSunrise == 1 & ismember(T_SOL.is_mast,mastCond{iMast}));
                    y = [T_SOL.SOL(s1Ids);T_SOL.SOL(s2Ids)];
                elseif iSubplot == 2
                    s1Ids = find(T_SOL.season == iSeason & T_SOL.isSunrise == 1 & ismember(T_SOL.is_mast,mastCond{iMast}));
                    s2Ids = find(T_SOL.season == kSeason & T_SOL.isSunrise == 1 & ismember(T_SOL.is_mast,mastCond{iMast}));
                    y = [mean([T_SOL.awakeIdx(s1Ids),T_SOL.asleepIdx(s1Ids)],2);mean([T_SOL.awakeIdx(s2Ids),T_SOL.asleepIdx(s2Ids)],2)];
                elseif iSubplot == 3
                    s1Ids = find(T_SOL.season == iSeason & T_SOL.isSunrise == 0 & ismember(T_SOL.is_mast,mastCond{iMast}));
                    s2Ids = find(T_SOL.season == kSeason & T_SOL.isSunrise == 0 & ismember(T_SOL.is_mast,mastCond{iMast}));
                    y = [T_SOL.SOL(s1Ids);T_SOL.SOL(s2Ids)];
                else
                    s1Ids = find(T_SOL.season == iSeason & T_SOL.isSunrise == 0 & ismember(T_SOL.is_mast,mastCond{iMast}));
                    s2Ids = find(T_SOL.season == kSeason & T_SOL.isSunrise == 0 & ismember(T_SOL.is_mast,mastCond{iMast}));
                    y = [mean([T_SOL.awakeIdx(s1Ids),T_SOL.asleepIdx(s1Ids)],2);mean([T_SOL.awakeIdx(s2Ids),T_SOL.asleepIdx(s2Ids)],2)];
                end
                if ismember(iMast,2:3)
                    mastSOL{iSubplot,iMast-1,iSeason} = y(1:numel(s1Ids));
                end
                group = [zeros(size(s1Ids));ones(size(s2Ids))];
                pMat_sol(iSeason,kSeason) = anova1(y,group,'off');
            end
        end
        disp(rowNames{iSubplot});
        flip(pMat_sol)
        writematrix(flip(pMat_sol),fullfile(exportPath,sprintf('%s_%s.csv',rowNames{iSubplot},mastNames{iMast})));
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
        set(gca,'fontsize',gcaFontSize);
        set(gca,'ydir','normal');
        drawnow;
    end
    cb = cbAside(gca,'p-value','k');
    cb.FontSize = gcaFontSize-2;
    cb.TickLabels = [0 0.05];
    
    if doSave
        %         print(gcf,'-painters','-depsc',fullfile(exportPath,[SOLFigureFiles{iMast},'.eps'])); % required for vector lines
        saveas(gcf,fullfile(exportPath,[SOLFigureFiles{iMast},'.jpg']),'jpg');
        close(gcf);
    end
end
%%
% do mast comparison
mastPmat = [];
for iSeason = 1:4
    for iCond = 1:4
        y = [mastSOL{iCond,1,iSeason};mastSOL{iCond,2,iSeason}];
        group = [zeros(size(mastSOL{iCond,1,iSeason}));ones(size(mastSOL{iCond,2,iSeason}))];
        mastPmat(iCond,iSeason) = anova1(y,group,'off');
    end
end
writematrix(mastPmat,fullfile(exportPath,'SOL_mastPmat.csv'));

%% !RUNFIG! plot all seasons with mast cond
mastCond = {[0,1],0,1};
mastTitle = {'',' - Non-mast',' - Mast'};
doMast = 0; % if 0, it assumes this is part of the big plot (all conds), else, creates the supp figure
ylims = [0.2,0.3,0.2,0.3,0.4];
if doMast
    close all;
end
for iMast = 1:3
    if doMast
        ff(500,900);
    else
        if iMast > 1
            continue;
        end
    end
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
        
        % LEFT COLUMN
        if doMast
            subplot(5,2,prc(2,[iSeason,1]));
        else
            subplot(rows,cols,figs{4+prc(2,[iSeason,1])});
            cla(gca);
        end
        binEdges = linspace(0,250,50);
        histogram(SOLs_sunrise,binEdges,'FaceColor',sunriseColor,'Normalization','probability');
        hold on;
        histogram(SOLs_sunset,binEdges,'FaceColor','k','Normalization','probability');
        set(gca,'fontsize',gcaFontSize);
        xlim([min(binEdges) max(binEdges)]);
        if iSeason == 5
            xlabel('Duration (minutes)');
        else
            xticklabels([]);
        end
        ylim([0 ylims(iSeason)]);
        yticks(ylim);
        ylabel('Probability','VerticalAlignment','top');
        if iSeason == 1
            title({'Latency',sprintf('%s%s',seasonLabel,mastTitle{iMast})});
        else
            title(sprintf('%s%s',seasonLabel,mastTitle{iMast}));
        end
        grid on;
        legend({'Sunrise','Sunset'},'fontsize',gcaFontSize-2);
        legend boxoff;
        
        % RIGHT COLUMN
        if doMast
            subplot(5,2,prc(2,[iSeason,2]));
        else
            subplot(rows,cols,figs{4+prc(2,[iSeason,2])});
            cla(gca);
        end
        binEdges = linspace(720-400,720+400,50);
        histogram(sunrise_means,binEdges,'FaceColor',sunriseColor,'Normalization','probability');
        hold on;
        histogram(sunset_means,binEdges,'FaceColor','k','Normalization','probability');
        set(gca,'fontsize',gcaFontSize);
        xticks(720-60*6:120:720+60*6);
        xticklabels({'-6','-4','-2','0','2','4','6'});
        if iSeason == 5
            xlabel('Rel. Time (hours)');
        else
            xticklabels([]);
        end
        ylim([0 ylims(iSeason)]);
        yticks(ylim);
        yticklabels([]);
        % %     ylabel('Frequency');
        if iSeason == 1
            title({'Onset/Offset',sprintf('%s%s',seasonLabel,mastTitle{iMast})});
        else
            title(sprintf('%s%s',seasonLabel,mastTitle{iMast}));
        end
        legend off;
        grid on;
        % %     legend({'QB-AB','AB-QB'},'fontsize',11,'autoupdate','off');
        % %     legend boxoff;
        xline(720,'k-');
        
        drawnow;
        hold off;
    end
    
    if doMast
        print(gcf,'-painters','-depsc',fullfile(exportPath,sprintf('%s-%s.eps','T_SOL_cleanHistograms',mastTitle{iMast}))); % required for vector lines
        saveas(gcf,fullfile(exportPath,sprintf('%s-%s.jpg','T_SOL_cleanHistograms',mastTitle{iMast})),'jpg');
        close(gcf);
    end
end
%% !RUNFIG!
print(gcf,'-painters','-depsc',fullfile(exportPath,'QBdur_SOL.eps')); % required for vector lines
saveas(gcf,fullfile(exportPath,'QBdur_SOL.jpg'),'jpg');
