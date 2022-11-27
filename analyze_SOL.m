% setup in predict_awake.m
% SOL_T table setup and visualization (1/3), see also: /Users/matt/Documents/MATLAB/KRSP/analyze_circCorrSleep.m
doSave = 1;
doDebug = 0;
gcaFontSize = 12;
colors = lines(5);

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
        for iSun = 1:2 % 1 = sunrise, 2 = sunset
            if doDebug
                subplot(rows,cols,iSun);
            end
            if iSun == 1
                shiftBy = 720-round(mean(secDay(Tss.sunrise(Tss.doy == sq_doys(iRec))),1)/60); % sunrise
            else
                shiftBy = 720-round(mean(secDay(Tss.sunset(Tss.doy == sq_doys(iRec))),1)/60); % sunset
            end
            shiftedAsleep = circshift(sq_asleep(iRec,:),shiftBy);
            shiftedNest = circshift(sq_nest(iRec,:),shiftBy);
            
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
                ylim([-0.1 1.1]);
                yticks([0,alpha,1-alpha,1]);
                xlabel('Time (min)');
                
                yyaxis right;
                ln4 = plot(shiftedNest,'color',[colors(5,:),0.4],'linewidth',2);
                ylabel('In Nest');
                ylim([-0.1 1.1]);
                yticks([0,1]);
                yticklabels({'Out','In'});
                set(gca,'ycolor',colors(5,:));
                yyaxis left;
                
                set(gca,'fontsize',gcaFontSize);
                title(sprintf('Rel. to %s, iRec = %04d',mnLabels{iSun},iRec));
            end
     
                legend([ln1,ln2,ln3,ln4],{'QB','normpdf','QB × normpdf','Nest'},'Autoupdate','off');
                m1 = [];
                m2 = [];
            end
            
            nestGradient = gradient(shiftedNest);
            nestExits = find(nestGradient == -0.5);
            nestEntrances = find(nestGradient == 0.5);
            if iSun == 1 % sunrise
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
                % nest exit
                [~,firstEE] = closest(nestExits,mean(asleepIdx,awakeIdx)); 
                if doDebug
                    plot([asleepIdx,awakeIdx],[asleepNorm(asleepIdx),asleepNorm(awakeIdx)],'color',repmat(0.5,[1,4]),'linewidth',10);
                    drawnow;
                end
            else % sunset
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
                % nest entrance
                [~,firstEE] = closest(nestEntrances,mean(asleepIdx,awakeIdx));
            end
            if doDebug
                %             plot([asleepIdx,awakeIdx],[0.5,0.5],'k-','linewidth',10);
                arrow([asleepIdx,0.5],[awakeIdx,0.5],'Length',5,'Ends',iSun);
                text(max([asleepIdx,awakeIdx])+10,0.5,sprintf('SOL = %i mins',abs(asleepIdx - awakeIdx)),...
                    'HorizontalAlignment','left','fontsize',gcaFontSize,'backgroundcolor','y');
                yyaxis right;
                plot(firstEE,iSun-1,'.','markerSize',40,'color',[colors(5,:),0.5]);
                plot(firstEE,iSun-1,'o','markerSize',20,'color',[colors(5,:),0.5],'linewidth',2);
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
            if isempty(firstEE)
                T_SOL.firstEE(iRow) = NaN;
            else
                T_SOL.firstEE(iRow) = firstEE;
            end
            T_SOL.doy(iRow) = sq_doys(iRec);
            for iSeason = 1:4
                if ismember(sq_doys(iRec),useDoys{iSeason})
                    T_SOL.season(iRow) = iSeason;
                end
            end
            T_SOL.is_mast(iRow) = ismember(sq_years(iRec),[2014,2019]);
            % % % %         T_SOL.asleepData(iRow) = {asleepNorm};
            T_SOL.sex(iRow) = sq_sex(iRec);
            T_SOL.is_preg(iRow) = sq_is_preg(iRec);
            T_SOL.sq_key_id(iRow) = sq_sqkeyrow(iRec);
        end
        if doSave && doDebug
            saveas(gcf,fullfile(debugPath,sprintf('iRec%04d.jpg',iRec)));
            close gcf;
        end
    end
    warning ('on','all');
    
    % rm outliers here, !!inspect for issues
    origSz = size(T_SOL,1);
    T_SOL(T_SOL.awakeIdx == 1,:) = [];
    T_SOL(T_SOL.asleepIdx == 1,:) = [];
    T_SOL(isnan(T_SOL.firstEE),:) = [];
    fnSize = size(T_SOL,1);
    fprintf('%i outliers removed (%1.1f%%), %i remain\n',origSz-fnSize,100*((origSz-fnSize)/origSz),fnSize);
    writetable(T_SOL,'T_SOL');
end
% use: T_SOL = readtable('T_SOL');

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
            pos = get(gca,'Position');
            set(gca,'Position',pos.*[1 1 1 0.97]);
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
            pos = get(gca,'Position');
            set(gca,'Position',pos.*[1 1 1 0.97]);
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
%% !RUNFIG! need to fix transparency manually
print(gcf,'-painters','-depsc',fullfile(exportPath,'QBdur_SOL.eps')); % required for vector lines
saveas(gcf,fullfile(exportPath,'QBdur_SOL.jpg'),'jpg');