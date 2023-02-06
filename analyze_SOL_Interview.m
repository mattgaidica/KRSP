colors = mycmap('/Users/mattgaidica/Documents/MATLAB/KRSP/util/seasons2.png',5);
useSex = 1;
mastCond = {[0,1]};
SOLs_sunrise = {};
SOLs_sunset = {};
sunrise_means = {};
sunset_means = {};
for iSeason = 1:5 % sseason 1 = all
    if iSeason == 1
        useIds = find(T_SOL.isSunrise==1 & ismember(T_SOL.sex,useSex));
        SOLs_sunrise{iSeason} = T_SOL.SOL(useIds);
        sunrise_means{iSeason} = mean([T_SOL.awakeIdx(useIds),T_SOL.asleepIdx(useIds)],2);
        
        useIds = find(T_SOL.isSunrise==0 & ismember(T_SOL.sex,useSex));
        SOLs_sunset{iSeason} = T_SOL.SOL(useIds);
        sunset_means{iSeason} = mean([T_SOL.awakeIdx(useIds),T_SOL.asleepIdx(useIds)],2);
    else
        useIds = find(T_SOL.isSunrise==1 & T_SOL.season == iSeason - 1 & ismember(T_SOL.sex,useSex));
        SOLs_sunrise{iSeason} = T_SOL.SOL(useIds);
        sunrise_means{iSeason} = mean([T_SOL.awakeIdx(useIds),T_SOL.asleepIdx(useIds)],2);
        
        useIds = find(T_SOL.isSunrise==0 & T_SOL.season == iSeason - 1 & ismember(T_SOL.sex,useSex));
        SOLs_sunset{iSeason} = T_SOL.SOL(useIds);
        sunset_means{iSeason} = mean([T_SOL.awakeIdx(useIds),T_SOL.asleepIdx(useIds)],2);
    end
end

%%
durationBinEdges = linspace(5,60,25);
timingBinEdges = linspace(720-400,720+400,50);
sunColors = [repmat(0.8,[1,3]);colors(1:4,:)];
seasonLabelsAll = [{'All'},seasonLabels(:)'];
sunLabels = {'Sunrise','Sunset'};
nestLabels = {'Exit','Entry'};

useSeasons = [2:5];
gcaFontSize = 16;
sunData = {sunrise_means,sunset_means};
close all;
ff(1000,400);
for iSS = 1:2
    subplot(1,2,iSS);
    lns = [];
    for iSeason = useSeasons
        lns(numel(lns)+1) = histogram(sunData{iSS}{iSeason},timingBinEdges,'FaceColor',sunColors(iSeason,:),...
            'Normalization','probability','FaceAlpha',0.5); %#ok<SAGROW> 
        hold on;
        set(gca,'fontsize',gcaFontSize);
        xticks(720-60*6:120:720+60*6);
        xticklabels({'-6','-4','-2','0','+2','+4','+6'});
        xlabel('Rel. Time (hours)');
        ylim([0 0.45]);
        yticks(ylim);
        ylabel('Probability');
        title(sprintf('Nest %s at %s',nestLabels{iSS},sunLabels{iSS}));
        grid on;
        xline(720,'k-');
        if numel(lns) == numel(useSeasons)
            legend(lns,seasonLabelsAll{useSeasons});
        end
    end
end
saveas(gcf,fullfile(exportPath,sprintf("sunriseSunsetDist_%s.jpg",strjoin(string(useSeasons),'-'))));

