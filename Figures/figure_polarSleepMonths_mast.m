%% mean
mastTitles = {'Mast','nMast'};
years_mast = [2014,2019]; % 2014,2019
years_nmast = [2015,2016,2017,2018,2020]; % 2015,2016,2017, *need 2018,2020
mast_years = {years_mast;years_nmast};

sIds = round(linspace(1,366,5));
seasonDoys = circshift(1:366,57);
useDoys = {seasonDoys(sIds(1):sIds(2)),seasonDoys(sIds(2):sIds(3)),...
    seasonDoys(sIds(3):sIds(4)),seasonDoys(sIds(4):sIds(5))};
monthNames = {'Winter','Spring','Summer','Autumn'};

close all
ff(1200,800);
iSubplot = 1;
% !! data should be circshift so sun is centered?
% or seasons should be year-centered, i.e. summer peaks at 171
for iSeason = 2:4
    mastIds = ismember(sq_doys,useDoys{iSeason}) & ismember(sq_years,years_mast);
    nmastIds = ismember(sq_doys,useDoys{iSeason}) & ismember(sq_years,years_nmast);

    monthData_nmast = imgaussfilt(mean(sq_asleep(nmastIds,:)),10,'padding','circular');
    monthData_mast = imgaussfilt(mean(sq_asleep(mastIds,:)),10,'padding','circular');

    subplot(2,3,prc(3,[1,iSubplot]));
    all_ps = [];
    for ii = 1:1440
        theseNmast = sq_asleep(nmastIds,ii);
        theseMast = sq_asleep(mastIds,ii);
        group = [zeros(size(theseNmast));ones(size(theseMast))];
        all_ps(ii) = anova1([theseNmast;theseMast],group,'off');
    end

    all_ps = pval_adjust(all_ps,'bonferroni');
    polarTime = linspace(-pi,pi,size(monthData_mast,2)) + pi/2;
    colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png',5);

    polarplot(polarTime,monthData_nmast,'color',colors(iSeason,:),'linewidth',7);
    hold on;
    % polarplot(polarTime,monthData_nmast + monthData_nmast_std,':','color','k','linewidth',1);
    polarplot(polarTime,monthData_mast,'color','k','linewidth',2);
    % polarplot(polarTime,monthData_nmast + monthData_mast_std,':','color','r','linewidth',1);

    pThresh = 0.001;
    polarplot(polarTime(all_ps<pThresh),ones(sum(all_ps<pThresh)),'r.','markersize',10);

    legend({'nmast','mast'},'location','southoutside');

    pax = gca;
    pax.ThetaTick = linspace(0,360,25);
    pax.ThetaTickLabel = 0:23;
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 14;
    pax.Layer = 'top';
    rlim([0 1]);
    rticks([]);
    pax.Color = [1 1 1];
    set(gca,'fontsize',14);

    title(sprintf("%s\nSleep Mean",monthNames{iSeason}));

    fs = 14;
    text(0,-1,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
    text(0.55,1.15,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));

    % var
    maxR = 1/3;
    monthData_nmast = imgaussfilt((maxR)-var(sq_asleep(nmastIds,:)),10,'padding','circular');
    monthData_mast = imgaussfilt((maxR)-var(sq_asleep(mastIds,:)),10,'padding','circular');

    subplot(2,3,prc(3,[2,iSubplot]));
    iSubplot = iSubplot + 1;
    all_ps = [];
    for ii = 1:1440
        theseNmast = sq_asleep(nmastIds,ii);
        theseMast = sq_asleep(mastIds,ii);
        group = [zeros(size(theseNmast));ones(size(theseMast))];
        all_ps(ii) = anova1([theseNmast;theseMast],group,'off');
    end

    all_ps = pval_adjust(all_ps,'bonferroni');
    polarTime = linspace(-pi,pi,size(monthData_mast,2)) + pi/2;

    polarplot(polarTime,monthData_nmast,'color',colors(iSeason,:),'linewidth',7);
    hold on;
    % polarplot(polarTime,monthData_nmast + monthData_nmast_std,':','color','k','linewidth',1);
    polarplot(polarTime,monthData_mast,'color','k','linewidth',2);
    % polarplot(polarTime,monthData_nmast + monthData_mast_std,':','color','r','linewidth',1);

    polarplot(polarTime(all_ps<pThresh),maxR*ones(sum(all_ps<pThresh)),'r.','markersize',10);

    legend({'nmast','mast'},'location','southoutside');

    pax = gca;
    pax.ThetaTick = linspace(0,360,25);
    pax.ThetaTickLabel = 0:23;
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 14;
    pax.Layer = 'top';
    rlim([0 maxR]);
    rticks([]);
    pax.Color = [1 1 1];
    set(gca,'fontsize',14);

    title(sprintf("%s\nSleep Consistency",monthNames{iSeason}));

    fs = 14;
    text(0,-maxR,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
    text(0.55,.38,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
end