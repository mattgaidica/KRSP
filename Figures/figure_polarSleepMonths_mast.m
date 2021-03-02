%% mean
% see also: /Users/matt/Documents/MATLAB/KRSP/run_this.m for comparing
% annual temp/light
weatherPath = '/Users/matt/Documents/Data/KRSP/HainesJunction_DailyTemps_Master.csv';
T_weather = readtable(weatherPath);

mastTitles = {'Mast','nMast'};
years_mast = [2014,2019]; % 2014,2019
years_nmast = [2015,2016,2017,2018,2020]; % 2015,2016,2017, *need 2018,2020
mast_years = {years_mast;years_nmast};

sIds = round(linspace(1,366,5));
seasonDoys = circshift(1:366,57); % centers at 171, so light is equal
% seasonDoys = circshift(1:366,57-21); % centers at 192, so temp is equal
useDoys = {seasonDoys(sIds(1):sIds(2)),seasonDoys(sIds(2):sIds(3)),...
    seasonDoys(sIds(3):sIds(4)),seasonDoys(sIds(4):sIds(5))};
monthNames = {'Winter','Spring','Summer','Autumn'};
linecolors = lines(5);

pThresh = 0.001;

% close all
ff(1200,800);
iSubplot = 1;
% data comes in sunrise-centered, this shifts it to solar-noon centered to
% eliminate bias of sunrise/sunset
for iSeason = 2:4
    mastIds = find(ismember(sq_doys,useDoys{iSeason}) & ismember(sq_years,years_mast) & ismember(sq_sex,1));
    nmastIds = find(ismember(sq_doys,useDoys{iSeason}) & ismember(sq_years,years_nmast) & ismember(sq_sex,1));
    
    % what if I only use days that appear in both mast and nmast years?
    newmastIds = [];
    newnmastIds = [];
    for iDoy = 1:366
        theseIdx_nmast = nmastIds(sq_doys(nmastIds)==iDoy);
        theseIdx_mast = mastIds(sq_doys(mastIds)==iDoy);
        onlyInclude = min([numel(theseIdx_mast),numel(theseIdx_nmast)]);
        if onlyInclude > 0
            if numel(theseIdx_mast) == numel(theseIdx_nmast)
                newmastIds = [newmastIds,theseIdx_mast];
                newnmastIds = [newnmastIds,theseIdx_nmast];
            elseif numel(theseIdx_mast) > numel(theseIdx_nmast)
                newmastIds = [newmastIds,theseIdx_mast(randi(numel(theseIdx_mast),[1,onlyInclude]))];
                newnmastIds = [newnmastIds,theseIdx_nmast];
            else
                newmastIds = [newmastIds,theseIdx_mast];
                newnmastIds = [newnmastIds,theseIdx_nmast(randi(numel(theseIdx_nmast),[1,onlyInclude]))];
            end
        end
    end
    
    mastIds = newmastIds;
    nmastIds = newnmastIds;
    
    % is the temperature the same these years?
    meanTemps_nmast = [];
    for ii = 1:numel(nmastIds)
        if sq_years(nmastIds(ii)) == 2020
            continue;
        end
        thisRow = find(T_weather.Year == sq_years(nmastIds(ii)) & T_weather.Julian_Date == sq_doys(nmastIds(ii)));
        if ~isnan(T_weather.Mean_Temp(thisRow))
            meanTemps_nmast = [meanTemps_nmast,T_weather.Mean_Temp(thisRow)];
        end
    end
    meanTemps_mast = [];
    for ii = 1:numel(mastIds)
        if sq_years(nmastIds(ii)) == 2020
            continue;
        end
        thisRow = find(T_weather.Year == sq_years(mastIds(ii)) & T_weather.Julian_Date == sq_doys(mastIds(ii)));
        if ~isnan(T_weather.Mean_Temp(thisRow))
            meanTemps_mast = [meanTemps_mast,T_weather.Mean_Temp(thisRow)];
        end
    end
    y = [meanTemps_nmast,meanTemps_mast];
    group = [zeros(size(meanTemps_nmast)),ones(size(meanTemps_mast))];
    p = anova1(y,group,'off');
    fprintf("season %i, temps p = %1.5f\n",iSeason,p);
    


    % this was sunsrise-centered
% % % %     monthData_nmast = imgaussfilt(mean(sq_asleep(nmastIds,:)),10,'padding','circular');
% % % %     monthData_mast = imgaussfilt(mean(sq_asleep(mastIds,:)),10,'padding','circular');
    
    nmastShifted = [];
    for iD = 1:numel(nmastIds)
        shiftBy = round(minutes(Tss.solar_noon(sq_doys(nmastIds(iD))) - Tss.sunrise(sq_doys(nmastIds(iD)))));
        nmastShifted(iD,:) = circshift(sq_asleep(nmastIds(iD),:),-shiftBy);
    end
    mastShifted = [];
    for iD = 1:numel(mastIds)
        shiftBy = round(minutes(Tss.solar_noon(sq_doys(mastIds(iD))) - Tss.sunrise(sq_doys(mastIds(iD)))));
        mastShifted(iD,:) = circshift(sq_asleep(mastIds(iD),:),-shiftBy);
    end
    
    % are sunrise times similarly distributed in this season?
    y = [secDay(Tss.sunrise(sq_doys(mastIds)));secDay(Tss.sunrise(sq_doys(nmastIds)))];
    group = [zeros(size(mastIds)),ones(size(nmastIds))];
    p = anova1(y,group,'off');
    fprintf("season %i, sunrise p = %1.5f\n",iSeason,p);
    
    morning = round(mean(secDay(Tss.solar_noon(useDoys{iSeason})))/60) - ...
        round(mean(secDay(Tss.sunrise(useDoys{iSeason})))/60); % minutes
    afternoon = round(mean(secDay(Tss.sunset(useDoys{iSeason})))/60) - ...
        round(mean(secDay(Tss.solar_noon(useDoys{iSeason})))/60); % minutes

    subplot(2,3,prc(3,[1,iSubplot]));
    all_ps = [];
    for ii = 1:1440
        theseNmast = nmastShifted(:,ii);
        theseMast = mastShifted(:,ii);
        group = [zeros(size(theseNmast));ones(size(theseMast))];
        all_ps(ii) = anova1([theseNmast;theseMast],group,'off');
    end
    all_ps = pval_adjust(all_ps,'bonferroni');
    
    sunBuffer = 0.1;
    polarTime = linspace(-pi,pi,1440);
    colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png',5);
    
    polarplot(polarTime,maxR*ones(1440,1),'color',repmat(0.7,[3,1]),'linewidth',0.75); % outer circle
    hold on;
    % actual data, colors(iSeason,:)
    lns = [];
    lns(1) = polarplot(polarTime,imgaussfilt(mean(nmastShifted),10,'padding','circular'),'color',repmat(0.5,[1,3]),'linewidth',7);
    lns(2) = polarplot(polarTime,imgaussfilt(mean(mastShifted),10,'padding','circular'),'color','k','linewidth',2);
    
    % mark morning/afternoon sun (and night)
    sunTime = [linspace(-(2*pi)*(morning/1440),0,morning) linspace(0,(2*pi)*(afternoon/1440),afternoon)];
    polarplot(polarTime,ones(1440,1)+sunBuffer,'color','k','linewidth',8); %night
    polarplot(sunTime,ones(numel(sunTime),1)+sunBuffer,'color',linecolors(3,:),'linewidth',8); %day

    polarplot(polarTime(all_ps<pThresh),ones(sum(all_ps<pThresh))+(sunBuffer/2),'r.','markersize',8);
    
    legend(lns,{'nmast','mast'},'location','southoutside');

    pax = gca;
    pax.ThetaTick = linspace(0,360,13);
    pax.ThetaTickLabel = {'Solar Noon','','','+6','','','±12','','','-6','',''};
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 14;
    pax.Layer = 'top';
    rlim([0 1+sunBuffer]);
    rticks([]);
    pax.Color = [1 1 1];
    set(gca,'fontsize',14);

    title(sprintf("%s\nSleep Mean",monthNames{iSeason}));

    fs = 14;
    text(0,1,["\uparrow","100%"],'horizontalalignment','center','verticalalignment','top','fontsize',fs,'color',repmat(0.4,[3,1]));

    
    % var
    maxR = 1/3;
    
    subplot(2,3,prc(3,[2,iSubplot]));
    iSubplot = iSubplot + 1;

    sunBuffer = 0.1 * maxR;
    polarplot(polarTime,ones(1440,1),'color',repmat(0.7,[3,1]),'linewidth',0.75); % outer circle
    hold on;
    % actual data
    lns = [];
    lns(1) = polarplot(polarTime,imgaussfilt(maxR-var(nmastShifted),10,'padding','circular'),'color',repmat(0.5,[1,3]),'linewidth',7);
    lns(2) = polarplot(polarTime,imgaussfilt(maxR-var(mastShifted),10,'padding','circular'),'color','k','linewidth',2);

    % mark morning/afternoon sun (and night)
    polarplot(polarTime,maxR*ones(1440,1)+sunBuffer,'color','k','linewidth',8); %night
    polarplot(sunTime,maxR*ones(numel(sunTime),1)+sunBuffer,'color',linecolors(3,:),'linewidth',8); %day

    polarplot(polarTime(all_ps<pThresh),maxR*ones(sum(all_ps<pThresh))+(sunBuffer/2),'r.','markersize',8);

    legend(lns,{'nmast','mast'},'location','southoutside');

    pax = gca;
    pax.ThetaTick = linspace(0,360,13);
    pax.ThetaTickLabel = {'Solar Noon','','','+6','','','±12','','','-6','',''};
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 14;
    pax.Layer = 'top';
    rlim([0 maxR+sunBuffer]);
    rticks([]);
    pax.Color = [1 1 1];
    set(gca,'fontsize',14);

    title(sprintf("%s\nSleep Consistency",monthNames{iSeason}));

    fs = 14;
    text(0,maxR,["\uparrow","100%"],'horizontalalignment','center','verticalalignment','top','fontsize',fs,'color',repmat(0.4,[3,1]));
end