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

iSeason = 3;
gcaFontSize = 16;
% close all;
% ff(1000,400);

subplot(121);
% histogram(SOLs_sunrise{iSeason},durationBinEdges,'FaceColor',sunColors(iSeason,:),'Normalization','probability');
counts = histcounts(SOLs_sunrise{iSeason},durationBinEdges,'Normalization','probability');
countsSm = smoothdata(equalVectors(counts,1000),'gaussian',100);
plot(countsSm,'linewidth',3,'color',sunColors(iSeason,:));
hold on;
hold on;
% histogram(SOLs_sunset{iSeason},durationBinEdges,'FaceColor','k','Normalization','probability');
set(gca,'fontsize',gcaFontSize);
% xlim([min(durationBinEdges) max(durationBinEdges)]);
xlabel('Duration (minutes)');
ylim([0 0.4]);
yticks(ylim);
ylabel('Probability','VerticalAlignment','top');
title(sprintf('%s',seasonLabelsAll{iSeason}));
grid on;
legend({'Sunrise','Sunset'},'fontsize',gcaFontSize-2);
legend boxoff;

subplot(122);
histogram(sunrise_means{iSeason},timingBinEdges,'FaceColor',sunColors(iSeason,:),'Normalization','probability');
hold on;
histogram(sunset_means{iSeason},timingBinEdges,'FaceColor','k','Normalization','probability');
set(gca,'fontsize',gcaFontSize);
xticks(720-60*6:120:720+60*6);
xticklabels({'-6','-4','-2','0','2','4','6'});
xlabel('Rel. Time (hours)');
ylim([0 0.4]);
yticks(ylim);
title(sprintf('%s%s',seasonLabelsAll{iSeason}));
legend off;
grid on;
legend({'Sunrise','Sunset'},'fontsize',gcaFontSize-2,'AutoUpdate','off');
legend boxoff;
xline(720,'k-');
