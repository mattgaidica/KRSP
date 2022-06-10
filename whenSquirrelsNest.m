% from: /Users/matt/Documents/MATLAB/KRSP/analyze_circCorrSleep.m
% question that Ben had...
useMonths =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
if ~exist('T_weather','var')
    weatherPath = '/Users/matt/Documents/Data/KRSP/HainesJunction_DailyTemps_Master.csv';
    T_weather = readtable(weatherPath);
end

close all
subplotMargins = [.1,0]; % [vert, horz]
doSave = true;
nHalfWindow = 30;
allDoys = 1:366;
op = 0.1;
ms = 20;
titleLabels = {'Hour of Day','Solar Noon'};
h = ff(900,500);
for iSun = 1:2
    subplot(1,2,iSun);
    sunBuffer = 0.05;
    sunrises = [];
    sunsets = [];
    theseAsleep_mean = [];
    for iDoy = 1:366
        shiftDoys = circshift(allDoys,-iDoy+1+nHalfWindow);
        theseDoys = shiftDoys(1:nHalfWindow*2+1);
        useIds = find(ismember(sq_doys,theseDoys));% & ismember(sq_sex,iSex));
        sunrises(iDoy) = round(mean(secDay(Tss.sunrise(ismember(Tss.doy,theseDoys))),1)/60); % minutes
        sunsets(iDoy) = round(mean(secDay(Tss.sunset(ismember(Tss.doy,theseDoys))),1)/60);
        if numel(useIds) > 1
            % data comes in uncentered
            if iSun == 1
                asleepDay_mean = mean(sq_nest(useIds,:));
            else
                asleepDay_mean = circshift(mean(sq_nest(useIds,:)),-round(mean(secDay(Tss.noon(iDoy)),1)/60));
            end
            theseAsleep_mean(iDoy,:) = imgaussfilt(asleepDay_mean,1,'padding','circular');
        else
            theseAsleep_mean(iDoy,:) = NaN(1,1440);
        end
    end

    colors = seasonColors(1:366);
    modR = 0;
    
    for iDoy = circshift(1:366,0) % for color overlay
        polarTime = linspace(0,2*pi,size(theseAsleep_mean,2));
        polarplot(polarTime,theseAsleep_mean(iDoy,:),'color',[colors(iDoy,:),op],'linewidth',2);
        hold on;
        ms = 15;

        if iDoy > 366/2
            if iDoy < (366/2) + 366/4
                modR = modR - .00027;
            else
                modR = modR + .00027;
            end
        else
            ms = 20;
        end
        if iSun == 1
            polarplot(polarTime(sunrises(iDoy)),1+modR+sunBuffer,'.','color',[colors(iDoy,:),op],'markersize',ms);
            polarplot(polarTime(sunsets(iDoy)),1+modR+sunBuffer,'.','color',[colors(iDoy,:),op],'markersize',ms);
        else
            sn = round(mean(secDay(Tss.noon(Tss.doy == iDoy)),1)/60);
            polarplot(polarTime(1440 - (sn - sunrises(iDoy))),1+modR+sunBuffer,'.','color',[colors(iDoy,:),op],'markersize',ms);
            polarplot(polarTime(sunsets(iDoy) - sn),1+modR+sunBuffer,'.','color',[colors(iDoy,:),op],'markersize',ms);
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
    c.TickLabels = useMonths;
    c.TickDirection = 'out';
    c.FontSize = 11;
    %     title(['Asleep Mean, ',titles{iSex+1}]);
    title(sprintf("In-nest Mean\n(%s)",titleLabels{iSun}));

    polarplot(polarTime,ones(size(polarTime)),':','color',repmat(0,[3,1]));
    fs = 14;
    if iSun == 1
        text(0,-1,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
    end
%     text(0.55,1.15,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
end

if doSave
    saveas(h,fullfile(exportPath,'circularNest'),'jpg');
    close(h);
end