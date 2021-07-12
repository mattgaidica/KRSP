% setup with /Users/matt/Documents/MATLAB/KRSP/predict_awake.m
useMonths =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
if ~exist('T_weather','var')
    weatherPath = '/Users/matt/Documents/Data/KRSP/HainesJunction_DailyTemps_Master.csv';
    T_weather = readtable(weatherPath);
end

%% simple correlation between amount of sleep and daylight
% the behavioral demands in spring/autumn reduce sleep
theseAsleep_mean = [];
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
            theseAsleep_mean(iDoy,:) = mean(sq_asleep(useIds,:));
            theseOdba(iDoy,:) = mean(sq_odba(useIds,:));
        else
            theseAsleep_mean(iDoy,:) = NaN(1,1440);
            theseOdba(iDoy,:) =  NaN(1,1440);
        end
        theseLight(iDoy) = mean(Tss.day_length(Tss.doy == iDoy),1);
        theseTemp(iDoy) = nanmean(T_weather.Mean_Temp(T_weather.Julian_Date == iDoy));
    end
    meanAsleep = inpaint_nans(mean(theseAsleep_mean,2),4);
    meanOdba = inpaint_nans(mean(theseOdba,2),4);
    plot(normalize(imgaussfilt(meanAsleep,10,'padding','circular')),'k','linewidth',2);
    hold on;
    plot(normalize(mean(theseAsleep_mean,2)),'k:','linewidth',0.5);
    plot(normalize(imgaussfilt(meanOdba,10,'padding','circular')),'g','linewidth',2);
    plot(normalize(mean(theseOdba,2)),'g:','linewidth',0.5);
    plot(normalize(theseLight),'y','linewidth',2);
    plot(normalize(theseTemp),'r','linewidth',2);
    title(titles{iSex+1});
    grid on;
    xlim([1 366]);
end

%% sleep mean (and old) consistency plots
close all
h = ff(800,400);
subplotMargins = [.15,0]; % [vert, horz]
doSave = true;
nHalfWindow = 30;
allDoys = 1:366;
op = 0.1;
ms = 20;
titleLabels = {'Hour of Day','Solar Noon'};

for iSun = 1:2
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
                asleepDay_mean = mean(sq_asleep(useIds,:));
            else
                asleepDay_mean = circshift(mean(sq_asleep(useIds,:)),-round(mean(secDay(Tss.noon(iDoy)),1)/60));
            end
            theseAsleep_mean(iDoy,:) = imgaussfilt(asleepDay_mean,1,'padding','circular');
        else
            theseAsleep_mean(iDoy,:) = NaN(1,1440);
        end
    end

    subplot(1,2,iSun);
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
    c.FontSize = 9;
    %     title(['Asleep Mean, ',titles{iSex+1}]);
    title({'Asleep Mean',titleLabels{iSun}});

    polarplot(polarTime,ones(size(polarTime)),':','color',repmat(0,[3,1]));
    fs = 14;
    if iSun == 1
        text(0,-1,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
    end
%     text(0.55,1.15,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));

    % consistency
% % % %     subplot(2,2,prc(2,[iSun,2]));
% % % %     maxR = 1/3;
% % % %     sunBuffer = sunBuffer * maxR;
% % % %     theseAsleep = [];
% % % %     for iDoy = 1:366
% % % %         shiftDoys = circshift(allDoys,-iDoy+1+nHalfWindow);
% % % %         useDoys = shiftDoys(1:nHalfWindow*2+1);
% % % %         useIds = find(ismember(sq_doys,useDoys));% & ismember(sq_sex,iSex));
% % % %         if numel(useIds) > 1
% % % %             if iSun == 1
% % % %                 monthData = maxR-var(sq_asleep(useIds,:));
% % % %             else
% % % %                 monthData = circshift(maxR-var(sq_asleep(useIds,:)),-round(secDay(Tss.solar_noon(iDoy))/60));
% % % %             end
% % % %             theseAsleep(iDoy,:) = imgaussfilt(monthData,1440/24,'padding','circular');
% % % %         else
% % % %             theseAsleep(iDoy,:) = NaN(1,1440);
% % % %         end
% % % %     end
% % % % 
% % % %     for iDoy = circshift(1:366,0) % for color overlay
% % % %         polarplot(polarTime,theseAsleep(iDoy,:),'color',[colors(iDoy,:),op],'linewidth',2);
% % % %         hold on;
% % % %         if iSun == 1
% % % %             ms = 15;
% % % %             if iDoy > 366/2
% % % %                 if iDoy < (366/2) + 366/4
% % % %                     modR = modR - (.00027*maxR);
% % % %                 else
% % % %                     modR = modR + (.00027*maxR);
% % % %                 end
% % % %             else
% % % %                 ms = 20;
% % % %             end
% % % %             polarplot(polarTime(round(sunrises(iDoy))),maxR+modR+sunBuffer,'.','color',[colors(iDoy,:),op],'markersize',ms);
% % % %             polarplot(polarTime(round(sunsets(iDoy))),maxR+modR+sunBuffer,'.','color',[colors(iDoy,:),op],'markersize',ms);
% % % %         end
% % % %     end
% % % %     pax = gca;
% % % %     pax.ThetaTick = linspace(0,360,25);
% % % %     pax.ThetaZeroLocation = 'top';
% % % %     pax.ThetaDir = 'clockwise';
% % % %     pax.FontSize = 14;
% % % %     pax.Layer = 'top';
% % % %     rlim([0 maxR+sunBuffer+(sunBuffer*.1)]);
% % % %     rticks([]);
% % % %     pax.Color = [1 1 1];
% % % %     if iSun == 1
% % % %         pax.ThetaTickLabel = 0:23;
% % % %     else
% % % %         pax.ThetaTickLabel = {'0','','','','','','+6','','','','','','±12','','','','','','-6','','','','','',};
% % % %     end
% % % % 
% % % %     c = colorbar('location','southoutside');
% % % %     colormap(colors);
% % % %     c.Limits = [0,1];
% % % %     c.Ticks = linspace(0,1,12);
% % % %     c.TickLabels = months;
% % % %     c.TickDirection = 'out';
% % % %     c.FontSize = 9;
% % % %     %     title(['Asleep Consistency, ',titles{iSex+1}]);
% % % %     title({'Asleep Consistency',titleLabels{iSun}});
% % % % 
% % % %     polarplot(polarTime,maxR*ones(size(polarTime)),':','color',repmat(0,[3,1]));
% % % %     fs = 14;
% % % %     if iSun == 1
% % % %         text(0,-maxR,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
% % % %     end
%     text(0.55,.38,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
end

% XY = [-.1,1;-.1,1];
% addFigureLabels(h,XY);
% setFig('','',2); % not sure if this is needed?
set(h,'PaperPositionMode','auto');

if doSave
%     print(gcf,'-painters','-depsc',fullfile(exportPath,'nestSleepOverlap.eps'));
    saveas(h,fullfile(exportPath,'circularSleep.eps'),'epsc');
    saveas(h,fullfile(exportPath,'circularSleep.jpg'),'jpg');
    close(h);
end

%% seasonal sleep w/ std
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);

close all
h = ff(800,900);
h_odba = ff(1400,900);
h_scat = ff(800,900);
t = linspace(0,24,1440);
for iSun = 1:2
    theseAsleep_mean = [];
    theseAsleep_std = [];
    
    theseOdba_mean = [];
    theseOdba_std = [];
    for iDoy = 1:366
        useIds = find(ismember(sq_doys,iDoy));% & ismember(sq_sex,iSex));
        if numel(useIds) > 1
            if iSun == 1
%                 monthData_mean = mean(sq_asleep(useIds,:));
%                 monthData_std = std(sq_asleep(useIds,:));
                asleepDay_mean = circshift(mean(sq_asleep(useIds,:),1),720-round(mean(secDay(Tss.sunrise(Tss.doy == iDoy)),1)/60));
                asleepDay_std = circshift(std(sq_asleep(useIds,:),[],1),720-round(mean(secDay(Tss.sunrise(Tss.doy == iDoy)),1)/60));
                
                odbaDay_mean = circshift(mean(sq_odba(useIds,:),1),720-round(mean(secDay(Tss.sunrise(Tss.doy == iDoy)),1)/60));
                odbaDay_std = circshift(std(sq_odba(useIds,:),[],1),720-round(mean(secDay(Tss.sunrise(Tss.doy == iDoy)),1)/60));
            else
                asleepDay_mean = circshift(mean(sq_asleep(useIds,:),1),720-round(mean(secDay(Tss.sunset(Tss.doy == iDoy)),1)/60));
                asleepDay_std = circshift(std(sq_asleep(useIds,:),[],1),720-round(mean(secDay(Tss.sunset(Tss.doy == iDoy)),1)/60));
                
                odbaDay_mean = circshift(mean(sq_odba(useIds,:),1),720-round(mean(secDay(Tss.sunset(Tss.doy == iDoy)),1)/60));
                odbaDay_std = circshift(std(sq_odba(useIds,:),[],1),720-round(mean(secDay(Tss.sunset(Tss.doy == iDoy)),1)/60));
            end
            theseAsleep_mean(iDoy,:) = asleepDay_mean;
            theseAsleep_std(iDoy,:) = asleepDay_std;
            
            theseOdba_mean(iDoy,:) = odbaDay_mean;
            theseOdba_std(iDoy,:) = odbaDay_std;
        else
            theseAsleep_mean(iDoy,:) = NaN(1,1440);
            theseAsleep_std(iDoy,:) = NaN(1,1440);
            
            theseOdba_mean(iDoy,:) = NaN(1,1440);
            theseOdba_std(iDoy,:) = NaN(1,1440);
        end
    end
    
    for iSeason = 1:4
        figure(h); % main mean fig
        subplot(4,2,prc(2,[iSeason,iSun]));
        x = nanmean(theseAsleep_mean(useDoys{iSeason},:));
        x_std = nanmean(theseAsleep_std(useDoys{iSeason},:));
        curve1 = x + x_std;
        curve2 = x - x_std;
        t2 = [t,fliplr(t)];
        fillArea = [curve1,fliplr(curve2)];
        fill(t2,fillArea,colors(iSeason,:),'FaceAlpha',0.15,'EdgeColor','none');
        hold on;
        plot(t,x,'color',colors(iSeason,:),'linewidth',3);
        
        xlim([0 24]);
        xticks(0:24);
        xticklabels({'±12','','','','','','-6','','','','','','0','','','','','','+6','','','','','±12'});
        xline(12,'k:');
        ylim([-0.25 1.25]);
        yticks([0 1]);
        ylabel('asleep mean')
        if iSun == 1
            xlabel('rel. to sunrise');
        else
            xlabel('rel. to sunset');
        end
        set(gca,'fontsize',14);
        title(seasonLabels{iSeason});
        drawnow;
        
        figure(h_odba);
        % mean
        subplot(4,4,prc(4,[iSeason,iSun*2-1]));
        plot(t,normalize(x,'range'),'color',colors(iSeason,:),'linewidth',2);
        hold on;
        xlim([0 24]);
        xticks(0:24);
        xticklabels({'±12','','','','','','-6','','','','','','0','','','','','','+6','','','','','±12'});
%         ylim([-0.25 1.25]);
        yticks([0 1]);
        ylabel('asleep mean')
        if iSun == 1
            xlabel('rel. to sunrise');
        else
            xlabel('rel. to sunset');
        end
        set(gca,'fontsize',14);
        title([seasonLabels{iSeason},' Mean']);
        
        x_odba = -nanmean(theseOdba_mean(useDoys{iSeason},:));
        plot(t,normalize(x_odba,'range'),'-','color',colors(iSeason,:).^2,'linewidth',2);
        if iSun == 1
            if iSeason == 1 || iSeason == 2
                legendLocation = 'southwest';
            else
                legendLocation = 'northeast';
            end
            legend({'asleep','odba^{-1}'},'location',legendLocation,'autoupdate','off');
            legend boxoff;
        end
        xline(12,'k:');
        
        % std
        subplot(4,4,prc(4,[iSeason,iSun*2]));
        plot(t,x_std - mean(x_std),'color',colors(iSeason,:),'linewidth',3);
        hold on;
        x_odba_std = -nanmean(theseOdba_std(useDoys{iSeason},:));
        thisStd = x_odba_std./(max(x_odba)-min(x_odba));
        plot(t,thisStd - mean(thisStd),'color',colors(iSeason,:).^2,'linewidth',3);
        xlim([0 24]);
        xticks(0:24);
        xticklabels({'±12','','','','','','-6','','','','','','0','','','','','','+6','','','','','±12'});
        xline(12,'k:');
        ylim([-1 1]);
        ylabel('std (a.u.)')
        if iSun == 1
            xlabel('rel. to sunrise');
        else
            xlabel('rel. to sunset');
        end
        set(gca,'fontsize',14);
        yticks([]);
        title([seasonLabels{iSeason},' Std']);
        drawnow;
        
        figure(h_scat);
        scatColors = [magma(720);flip(magma(720))];
        if iSun == 1
            subplot(2,2,iSeason);
            [r,p] = corr(normalize(x_odba,'range')',normalize(x,'range')');
            scatter(normalize(x_odba,'range'),normalize(x,'range'),10,scatColors,'filled');
            xlabel('odba^{-1} (norm.)');ylabel('asleep (norm.)');
            title({[seasonLabels{iSeason},' - Sunrise'],...
                sprintf('r = %1.2f, p = %1.2e',r,p)});
            c = colorbar('location','southoutside');
            c.Ticks = [0 0.5 1];
            c.TickLabels = {'±12','sunrise','±12'};
            colormap(scatColors);
            xlim([0 1]);
            ylim(xlim);
            xticks(0:0.2:1);
            yticks(xticks);
            set(gca,'fontsize',14);
        end
    end
end