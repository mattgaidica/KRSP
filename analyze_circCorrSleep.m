% setup with /Users/matt/Documents/MATLAB/KRSP/Figures/figure_rhythmicSleep_noDiff.m
months =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
if ~exist('T_weather','var')
    weatherPath = '/Users/matt/Documents/Data/KRSP/HainesJunction_DailyTemps_Master.csv';
    T_weather = readtable(weatherPath);
    ssPath = '/Users/matt/Documents/Data/KRSP/SunriseSunset';
    files = dir(fullfile(ssPath,'*.txt'));
    T_ss = readtable(fullfile(ssPath,files(4).name)); % 366 day year (simplify for now)
end
% temps(iSq) = T_weather.Mean_Temp(T_weather.Year == 2016 & T_weather.Julian_Date == thisDoy);

%% simple correlation between amount of sleep and daylight
% the behavioral demands in spring/autumn reduce sleep
theseAsleep = [];
theseLight = [];
theseTemp = [];
theseOdba = [];
close all
ff(1000,800);
titles = {'Female','Male'};
for iSex = 0:1
    subplot(2,1,iSex+1);
    for iDoy = 1:366
        useIds = find(ismember(sq_doys,iDoy) & ismember(sq_sex,iSex));
        if numel(useIds) > 1
            theseAsleep(iDoy,:) = mean(sq_asleep(useIds,:));
            theseOdba(iDoy,:) = mean(sq_odba(useIds,:));
        else
            theseAsleep(iDoy,:) = NaN(1,1440);
            theseOdba(iDoy,:) =  NaN(1,1440);
        end
        theseLight(iDoy) = T_ss.day_length(iDoy);
        theseTemp(iDoy) = nanmean(T_weather.Mean_Temp(T_weather.Julian_Date == iDoy));
    end
    meanAsleep = inpaint_nans(mean(theseAsleep,2),4);
    meanOdba = inpaint_nans(mean(theseOdba,2),4);
    plot(normalize(imgaussfilt(meanAsleep,10,'padding','circular')),'k','linewidth',2);
    hold on;
    plot(normalize(mean(theseAsleep,2)),'k:','linewidth',0.5);
    plot(normalize(imgaussfilt(meanOdba,10,'padding','circular')),'g','linewidth',2);
    plot(normalize(mean(theseOdba,2)),'g:','linewidth',0.5);
    plot(normalize(theseLight),'y','linewidth',2);
    plot(normalize(theseTemp),'r','linewidth',2);
    title(titles{iSex+1});
    grid on;
    xlim([1 366]);
end

%% investigate nap phenotype
% sleep mean by SEX
close all
ff(800,800);
subplot(1,2,1);
nHalfWindow = 30;
allDoys = 1:366;
op = 0.1;
ms = 20;
titleLabels = {'Hour of Day','Solar Noon'};

iSubplot = 1;
for iSun = 1:2
    sunBuffer = 0.05;
    sunrises = [];
    sunsets = [];
    theseAsleep = [];
    for iDoy = 1:366
        shiftDoys = circshift(allDoys,-iDoy+1+nHalfWindow);
        useDoys = shiftDoys(1:nHalfWindow*2+1);
        useIds = find(ismember(sq_doys,useDoys));% & ismember(sq_sex,iSex));
        sunrises(iDoy) = round(mean(secDay(T_ss.sunrise(useDoys))/60)); % minutes
        sunsets(iDoy) = round(mean(secDay(T_ss.sunset(useDoys))/60));
        if numel(useIds) > 1
            % data comes in sunrise centered
            if iSun == 1
                monthData = mean(sq_asleep(useIds,:));
            else
                monthData = circshift(mean(sq_asleep(useIds,:)),-round(secDay(T_ss.solar_noon(iDoy))/60));
            end
            theseAsleep(iDoy,:) = imgaussfilt(monthData,1440/24,'padding','circular');
        else
            theseAsleep(iDoy,:) = NaN(1,1440);
        end
    end
    
    if false
        ff(1000,500);
        for iDoy = 1:366
            plot(circshift(theseAsleep(iDoy,:),720),'color',[colors(iDoy,:),op],'linewidth',2);
            hold on;
        end
        xlim([1 1440]);
    end

    subplot(2,2,prc(2,[iSun,1]));
    polarTime = linspace(0,2*pi,size(theseAsleep,2));
    colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png',366);
    modR = 0;
    for iDoy = circshift(1:366,0) % for color overlay
        polarplot(polarTime,theseAsleep(iDoy,:),'color',[colors(iDoy,:),op],'linewidth',2);
        hold on;
        ms = 15;
        if iSun == 1
            if iDoy > 366/2
                if iDoy < (366/2) + 366/4
                    modR = modR - .00027;
                else
                    modR = modR + .00027;
                end
            else
                ms = 20;
            end
            polarplot(polarTime(sunrises(iDoy)),1+modR+sunBuffer,'.','color',[colors(iDoy,:),op],'markersize',ms);
            polarplot(polarTime(sunsets(iDoy)),1+modR+sunBuffer,'.','color',[colors(iDoy,:),op],'markersize',ms);
        end
    end
    pax = gca;
    pax.ThetaTick = linspace(0,360,25);
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 14;
    pax.Layer = 'top';
    rlim([0 1+sunBuffer+(sunBuffer*.1)]);
    rticks([]);
    pax.Color = [1 1 1];
    if iSun == 1
        pax.ThetaTickLabel = 0:23;
    else
        pax.ThetaTickLabel = {'0','','','','','','+6','','','','','','±12','','','','','','-6','','','','','',};
    end

    c = colorbar('location','southoutside');
    colormap(colors);
    c.Limits = [0,1];
    c.Ticks = linspace(0,1,12);
    c.TickLabels = months;
    c.TickDirection = 'out';
    c.FontSize = 9;
    %     title(['Asleep Mean, ',titles{iSex+1}]);
    title({'Asleep Mean',titleLabels{iSun}});

    polarplot(polarTime,ones(size(polarTime)),':','color',repmat(0,[3,1]));
    fs = 14;
    if iSun == 1
        text(0,-1,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
    end
%     text(0.55,1.15,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));

    subplot(2,2,prc(2,[iSun,2]));
    maxR = 1/3;
    sunBuffer = sunBuffer * maxR;

    % consistency
    theseAsleep = [];
    for iDoy = 1:366
        shiftDoys = circshift(allDoys,-iDoy+1+nHalfWindow);
        useDoys = shiftDoys(1:nHalfWindow*2+1);
        useIds = find(ismember(sq_doys,useDoys));% & ismember(sq_sex,iSex));
        if numel(useIds) > 1
            if iSun == 1
                monthData = maxR-var(sq_asleep(useIds,:));
            else
                monthData = circshift(maxR-var(sq_asleep(useIds,:)),-round(secDay(T_ss.solar_noon(iDoy))/60));
            end
            theseAsleep(iDoy,:) = imgaussfilt(monthData,1440/24,'padding','circular');
        else
            theseAsleep(iDoy,:) = NaN(1,1440);
        end
    end

    for iDoy = circshift(1:366,0) % for color overlay
        polarplot(polarTime,theseAsleep(iDoy,:),'color',[colors(iDoy,:),op],'linewidth',2);
        hold on;
        if iSun == 1
            ms = 15;
            if iDoy > 366/2
                if iDoy < (366/2) + 366/4
                    modR = modR - (.00027*maxR);
                else
                    modR = modR + (.00027*maxR);
                end
            else
                ms = 20;
            end
            polarplot(polarTime(round(sunrises(iDoy))),maxR+modR+sunBuffer,'.','color',[colors(iDoy,:),op],'markersize',ms);
            polarplot(polarTime(round(sunsets(iDoy))),maxR+modR+sunBuffer,'.','color',[colors(iDoy,:),op],'markersize',ms);
        end
    end
    pax = gca;
    pax.ThetaTick = linspace(0,360,25);
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 14;
    pax.Layer = 'top';
    rlim([0 maxR+sunBuffer+(sunBuffer*.1)]);
    rticks([]);
    pax.Color = [1 1 1];
    if iSun == 1
        pax.ThetaTickLabel = 0:23;
    else
        pax.ThetaTickLabel = {'0','','','','','','+6','','','','','','±12','','','','','','-6','','','','','',};
    end

    c = colorbar('location','southoutside');
    colormap(colors);
    c.Limits = [0,1];
    c.Ticks = linspace(0,1,12);
    c.TickLabels = months;
    c.TickDirection = 'out';
    c.FontSize = 9;
    %     title(['Asleep Consistency, ',titles{iSex+1}]);
    title({'Asleep Consistency',titleLabels{iSun}});

    polarplot(polarTime,maxR*ones(size(polarTime)),':','color',repmat(0,[3,1]));
    fs = 14;
    if iSun == 1
        text(0,-maxR,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
    end
%     text(0.55,.38,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
end