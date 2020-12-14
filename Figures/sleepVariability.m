% init with /Users/matt/Documents/MATLAB/KRSP/Figures/figure_rhythmicSleep_noDiff.m
% ^ stop after filtIds is set
tSeasons = round(linspace(1,366,5));
seasonDoys = circshift(1:366,60);
unSqs = unique(sq_ids);
titleLabels = {'winter','spring','summer','fall'};

Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
Tss_doys = day(Tss.sunrise,'dayofyear');
sunsetOffset = [];
for iDoy = 1:366
    afterSunrise = Tss.day_length(iDoy)/60; % minutes
    
    if afterSunrise > 720
        sunsetOffset(iDoy) = afterSunrise - 1440;
    else
        sunsetOffset(iDoy) = afterSunrise;
    end
end

odbaArr = NaN(366,1440);
odbaStdArr = NaN(366,1440);
asleepArr = NaN(366,1440);
asleepStdArr = NaN(366,1440);
for iDoy = 1:366
    useIds = ismember(sq_doys,iDoy) & filtIds;
    if sum(useIds) < 7
        continue;
    end
    odbaArr(iDoy,:) = mean(sq_odba_z(useIds,:));
    odbaStdArr(iDoy,:) = std(sq_odba_z(useIds,:));
    asleepArr(iDoy,:) = mean(sq_asleep(useIds,:));
    asleepStdArr(iDoy,:) = var(sq_asleep(useIds,:)) ./ sqrt(sum(useIds));
end

%% not used
nG = 1;
close all
ff(800,350);

subplot(221);
imagesc(imgaussfilt(odbaArr',nG));
colormap(magma);
title('odba mean');
caxis([0 3]);
set(gca,'ydir','normal');
colorbar;

subplot(223);
imagesc(imgaussfilt(odbaStdArr',nG));
colormap(magma);
title('odba std');
caxis([0 5]);
set(gca,'ydir','normal');
colorbar;

subplot(222);
imagesc(imgaussfilt(asleepArr',nG));
colormap(magma);
title('asleep mean');
caxis([0 1]);
set(gca,'ydir','normal');
colorbar;

subplot(224);
imagesc(imgaussfilt(asleepStdArr',nG));
colormap(magma);
title('asleep std');
caxis([0.2 0.5]);
set(gca,'ydir','normal');
colorbar;

%% main heat plots
lw = 2;
t = linspace(-720,720,size(asleepStdArr,2));
nG = 0.5;
op = 0.7;
sunsetOffset_nan = sunsetOffset;
sunsetOffset_nan(abs(diff(sunsetOffset)) > 50) = NaN;
% close all
ff(500,800,2);

subplot(411);
imagesc(1:366,t,imgaussfilt(odbaArr',nG),'AlphaData',~isnan(imgaussfilt(odbaArr',nG)));
colormap(magma);
caxis([-0.25 3]);
set(gca,'ydir','normal');
colorbar;
set(gca,'fontsize',14);
yticks([-720 -360 0 360 720]);
ylabel('Z_{time} (min)');
hold on;
plot(sunsetOffset_nan,'-','linewidth',lw,'color',[0 0 0 op]);
plot(xlim,[0 0],'--','linewidth',lw,'color',[0 0 0 op]);
xticks(round(linspace(1,366,5)));
title('Mean zODBA');

subplot(412);
imagesc(1:366,t,imgaussfilt(asleepArr',nG),'AlphaData',~isnan(imgaussfilt(asleepArr',nG)));
colormap(magma);
% caxis([0 5]);
set(gca,'ydir','normal');
colorbar;
set(gca,'fontsize',14);
yticks([-720 -360 0 360 720]);
ylabel('Z_{time} (min)');
hold on;
plot(sunsetOffset_nan,'-','linewidth',lw,'color',[0 0 0 op]);
plot(xlim,[0 0],'--','linewidth',lw,'color',[0 0 0 op]);
xticks(round(linspace(1,366,5)));
title('Mean Asleep');

subplot(413);
imagesc(1:366,t,imgaussfilt(asleepStdArr.^2',nG),'AlphaData',~isnan(imgaussfilt(asleepStdArr.^2',nG)));
colormap(magma);
title('Asleep Variability');
caxis([0 .07]);
set(gca,'ydir','normal');
colorbar;
set(gca,'fontsize',14);
yticks([-720 -360 0 360 720]);
ylabel('Z_{time} (min)');
hold on;
plot(sunsetOffset_nan,'-','linewidth',lw,'color',[0 0 0 op]);
plot(xlim,[0 0],'--','linewidth',lw,'color',[0 0 0 op]);
xticks(round(linspace(1,366,5)));

subplot(414);
asleepVar = asleepStdArr.^2;
% asleepFill = inpaint_nans(asleepVar,5);
asleepFill = fillmissing(asleepVar,'movmean',50);
filtStd = imgaussfilt(asleepFill,12,'padding','circular'); % variance
contourf(1:366,t,filtStd',5);
colormap(magma);
caxis([0 0.05]);
colorbar;
set(gca,'fontsize',14);
xlabel('day of year');
hold on;
plot(sunsetOffset_nan,'-','linewidth',lw,'color',[0 0 0 op]);
plot(xlim,[0 0],'--','linewidth',lw,'color',[0 0 0 op]);
yticks([-720 -360 0 360 720]);
ylabel('Z_{time} (min)');
title('Asleep Variability Topography');
xticks(round(linspace(1,366,5)));

%% circular plots
monthNames = {'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
tSpan = 30;
titleLabels = {'-30 min','Sunrise','+30 min';'-30 min','Sunset','+30 min'};
nS = 5;
polarVar_fill = imgaussfilt(asleepFill,nS,'padding','circular');
polarVar = polarVar_fill;
polarVar(isnan(asleepVar)) = NaN;
% night versions, just shift
dl = round(Tss.day_length/60);
asleepVar_n = asleepVar;
for ii = 1:366
    asleepVar_n(ii,:) = circshift(asleepVar_n(ii,:),-dl(ii));
end
asleepFill_n = fillmissing(asleepVar_n,'movmean',50); % make sure matches top
polarVar_n_fill = imgaussfilt(asleepFill_n,nS,'padding','circular'); % variance
polarVar_n = polarVar_n_fill; % use smoothed version
polarVar_n(isnan(asleepVar_n)) = NaN; % but put back NaN

ts = [720-tSpan,720,720+tSpan];
theta = linspace(0,2*pi,366);
close all
ff(900,600,2);
for iPlot = 1:3
    subplot(2,3,iPlot);
    v = sort(mean(polarVar_fill));
    loP = v(floor(366*.05));
    hiP = v(end-ceil(366*.05));
    polarplot(theta,repmat(loP,[1,366]),'color','r','linewidth',0.75,'linestyle',':');
    hold on;
    polarplot(theta,repmat(hiP,[1,366]),'color','r','linewidth',0.75,'linestyle',':');
    polarplot(theta,polarVar_fill(:,ts(iPlot)),'color',repmat(0.5,[1,3]),'linewidth',1);
    polarplot(theta,polarVar(:,ts(iPlot)),'color','k','linewidth',3);
    pax = gca;
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 14;
    pax.Layer = 'top';
    rlim([0 0.08]);
    rticks([]);
    pax.Color = [1 1 1];
    thetaticklabels(monthNames);
    title(titleLabels{1,iPlot});
    
    subplot(2,3,iPlot+3);
    v = sort(mean(polarVar_n_fill));
    loP = v(floor(366*.05));
    hiP = v(end-ceil(366*.05));
    polarplot(theta,repmat(loP,[1,366]),'color','r','linewidth',0.75,'linestyle',':');
    hold on;
    polarplot(theta,repmat(hiP,[1,366]),'color','r','linewidth',0.75,'linestyle',':');
    polarplot(theta,polarVar_n_fill(:,ts(iPlot)),'color',repmat(0.5,[1,3]),'linewidth',1);
    polarplot(theta,polarVar_n(:,ts(iPlot)),'color','k','linewidth',3);
    pax = gca;
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 14;
    pax.Layer = 'top';
    rlim([0 0.08]);
    rticks([]);
    pax.Color = [1 1 1];
    thetaticklabels(monthNames);
    title(titleLabels{2,iPlot});
end

%% inspect means
% close all
ff(400,900);
for iSeason = 1:4
    subplot(4,1,iSeason);
    disp(iSeason);
    useIds = ismember(sq_doys,seasonDoys(tSeasons(iSeason):tSeasons(iSeason+1)));
    plot(mean(sq_asleep(useIds,:)));
    ylim([0 1]);
    yyaxis right;
    plot(std(sq_asleep(useIds,:)));
    title(titleLabels{iSeason});
    ylim([0.35 0.55]);
    xlim([1 size(sq_odba_z,2)]);
end

%%
% close all
ff(400,900);
for iSeason = 1:4
    subplot(4,1,iSeason);
    disp(iSeason);
    useIds = ismember(sq_doys,seasonDoys(tSeasons(iSeason):tSeasons(iSeason+1)));
    plot(mean(sq_odba_z(useIds,:)));
    ylim([-0.25 4]);
    yyaxis right;
    plot(std(sq_odba_z(useIds,:)));
    title(titleLabels{iSeason});
    ylim([0 5]);
    xlim([1 size(sq_odba_z,2)]);
end