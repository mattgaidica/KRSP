filespath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
nBins = 96;
binEdges = linspace(0,86400-60,nBins+1);
doFig = false;
if do
%     files = dir(fullfile(filespath,'*.mat'));
    Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
    Tss_doys = day(Tss.sunrise,'dayofyear');
    all_hists = cell(366,1);
    close all
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq})
            load(fullfile(filespath,sqkey.filename{iSq})); % T, Tstat
% % % %             if ~isValidT(T,false)
% % % %                 disp(['Skipping ',sqkey.filename{iSq}]);
% % % %                 continue;
% % % %             end
        end
        T = detect_sleepWake(T,2);
        doys_all = day(T.datetime,'dayofyear');
        doys_unique = unique(doys_all);
        histArr = [];
        histNorm = [];
        for iDoy = 1:numel(doys_unique)
            useIds = find(doys_all == doys_unique(iDoy));
            useTimes = T.datetime(useIds);
            awakeArr = secDay(useTimes(T.awake(useIds) == 1)); % awake
            histArr(iDoy,:) = histcounts(awakeArr,binEdges);
            histNorm(iDoy,:) = histcounts(secDay(T.datetime(useIds(1))):60:...
                secDay(T.datetime(useIds(end))),binEdges);
            thisMat = all_hists{doys_unique(iDoy)};
            if isempty(thisMat)
                all_hists{doys_unique(iDoy)} = histArr(iDoy,:) ./ histNorm(iDoy,:);
            else
                thisMat(size(thisMat,1)+1,:) = histArr(iDoy,:) ./ histNorm(iDoy,:);
                all_hists{doys_unique(iDoy)} = thisMat;
            end
        end
        if doFig
            sunrise = secDay(Tss.sunrise(Tss_doys == doys_unique(iDoy)));
            sunset = secDay(Tss.sunset(Tss_doys == doys_unique(iDoy)));
            h = ff(1200,600);
            bar(linspace(0,86400,size(histArr,2)),sum(histArr,1)./sum(histNorm,1),'k');
            hold on;
            colors = lines(3);
            op = 0.5;
            lw = 8;
            plot([sunrise,sunrise],ylim,'-','linewidth',lw,'color',[colors(3,:),op]);
            plot([sunset,sunset],ylim,'-','linewidth',lw,'color',[0,0,0,op]);
            xticklabels(compose('%1.1f',24*xticks/(86400-60)));
            xlabel('hour of day');
            ylabel('frac. awake');
            title(sprintf('%s - %i days',sqkey.filename{iSq},size(histArr,1)),'interpreter','none');
            set(gca,'fontsize',14);
            saveas(h,strrep(fullfile(filespath,'_figs',sqkey.filename{iSq}),'.mat','.jpeg'));
            close(h);
        end
    end
    do = false;
end

mean_arr = [];
for iDoy = 1:numel(all_hists)
    mean_arr(iDoy,:) = mean(all_hists{iDoy});
end
close all
colors = lines(3);
ff(1200,600);
imagesc(mean_arr');
set(gca,'fontsize',14);
xlabel('day of year');
set(gca,'ydir','normal');
yticks(linspace(1,nBins,12));
yticklabels(round(linspace(0,24,numel(yticks))));
ylabel('hour of day');
c = colorbar;
caxis([0 1]);
c.Label.String = 'p(awake)';
c.Label.FontSize = 14;
colormap(magma);
hold on;
plot(secDay(Tss.sunrise)/(86400/nBins),'-','linewidth',4,'color',[colors(3,:),op]);
plot(secDay(Tss.sunset)/(86400/nBins),'-','linewidth',4,'color',[0,0,0,op]);