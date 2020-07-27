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
    asleepStdArr(iDoy,:) = var(sq_asleep(useIds,:));
end

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

%%
t = linspace(-720,720,size(asleepStdArr,2));
nG = 0.5;
op = 0.5;
sunsetOffset_nan = sunsetOffset;
sunsetOffset_nan(abs(diff(sunsetOffset)) > 50) = NaN;
close all
ff(500,800,2);

subplot(411);
imagesc(1:366,t,imgaussfilt(odbaArr',nG));
colormap(magma);
caxis([0 3]);
set(gca,'ydir','normal');
% colorbar;
set(gca,'fontsize',14);
yticks([-720 -360 0 360 720]);
ylabel('Z_{time} (min)');
hold on;
plot(sunsetOffset_nan,'-','linewidth',lw,'color',[0 0 0 op]);
plot(xlim,[0 0],'--','linewidth',lw,'color',[0 0 0 op]);
xticks(round(linspace(1,366,5)));
title('Mean ODBA');

subplot(412);
imagesc(1:366,t,imgaussfilt(asleepArr',nG));
colormap(magma);
% caxis([0 5]);
set(gca,'ydir','normal');
% colorbar;
set(gca,'fontsize',14);
yticks([-720 -360 0 360 720]);
ylabel('Z_{time} (min)');
hold on;
plot(sunsetOffset_nan,'-','linewidth',lw,'color',[0 0 0 op]);
plot(xlim,[0 0],'--','linewidth',lw,'color',[0 0 0 op]);
xticks(round(linspace(1,366,5)));
title('Mean Asleep');

subplot(413);
imagesc(1:366,t,imgaussfilt(asleepStdArr.^2',nG));
colormap(magma);
title('Asleep Variability');
% caxis([0.2 0.5]);
set(gca,'ydir','normal');
% colorbar;
set(gca,'fontsize',14);
yticks([-720 -360 0 360 720]);
ylabel('Z_{time} (min)');
hold on;
plot(sunsetOffset_nan,'-','linewidth',lw,'color',[0 0 0 op]);
plot(xlim,[0 0],'--','linewidth',lw,'color',[0 0 0 op]);
xticks(round(linspace(1,366,5)));

subplot(414);
filtStd = imgaussfilt(asleepStdArr.^2',8); % variance
filtStd_fill = inpaint_nans(filtStd,1);
nLevels = 5;
contourf(1:366,t,filtStd_fill,5);
colormap(magma);
% caxis([0.19 0.43]);
% colorbar;
set(gca,'fontsize',14);
xlabel('day of year');
hold on;
plot(sunsetOffset_nan,'-','linewidth',lw,'color',[0 0 0 op]);
plot(xlim,[0 0],'--','linewidth',lw,'color',[0 0 0 op]);
yticks([-720 -360 0 360 720]);
ylabel('Z_{time} (min)');
title('Asleep Variability Topography');
xticks(round(linspace(1,366,5)));
%%
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