% close all
ff(1200,800);
titleLabels = {'Sunrise SOL','Sunset SOL'};

reqSessions = 3;
for iSun = 1:2
    subplot(2,1,iSun);
    SOL_med = NaN(1,366);
    SOL_std = NaN(1,366);
    windowHalfSize = 1; % days
    for iDoy = 1:366
        doyArr = circshift(1:366,-iDoy+windowHalfSize);
        useDoys = doyArr(1:windowHalfSize*2);
        useIds = find(T_SOL.isSunrise == iSun-1 & ismember(T_SOL.doy,useDoys));
        if numel(useIds) > reqSessions - 1
            SOL_med(iDoy) = median(T_SOL.SOL(useIds),1);
            SOL_std(iDoy) = std(T_SOL.SOL(useIds),1);
        end
    end
    plot(SOL_med,'k','linewidth',2);
    hold on;
    plot(SOL_std,'r:');
    xlabel('doy');
    ylabel('minutes');
    title(titleLabels{iSun});
    xlim([1,366]);
    grid on;
end

%%
close all
useMonths =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
titleLabels_sol = {'Sunrise SOL','Sunset SOL'};
titleLabels_sleepWake = {'Sleep/Wake Timing'};
ff(1000,400);
reqSessions = 3;
rows = 1;
cols = 3;
colors = lines(7);
sleepWakeColors = colors([5,7],:);
seasonColors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
seasonCell = {seasonDoys(sIds(1):sIds(2)),seasonDoys(sIds(2):sIds(3)),...
    seasonDoys(sIds(3):sIds(4)),seasonDoys(sIds(4):sIds(5))};
sunriseKey = [1,0]; % goes: sunrise, sunset
for iSun = 1:2
    subplot(rows,cols,prc(cols,[1,iSun]));
    SOL_med = NaN(1,366);
    SOL_std = NaN(1,366);
    sleepWake_med = NaN(1,366);
    windowHalfSize = 1; % days
    for iDoy = 1:366
        doyArr = circshift(1:366,-iDoy+windowHalfSize);
        useDoys = doyArr(1:windowHalfSize*2);
        useIds = find(T_SOL.isSunrise == sunriseKey(iSun) & ismember(T_SOL.doy,useDoys));
        if numel(useIds) > reqSessions - 1
            SOL_med(iDoy) = median(T_SOL.SOL(useIds),1);
            SOL_std(iDoy) = std(T_SOL.SOL(useIds),1);
            sleepWake_med(iDoy) = median(mean([T_SOL.awakeIdx(useIds),T_SOL.asleepIdx(useIds)],2));
        end
    end
    theta = linspace(0,2*pi,366);
    %     polarplot(theta,SOL_med,'color',repmat(0.7,[4,1]),'linewidth',1);
    %     hold on;
    smoothMedSOL = imgaussfilt(SOL_med,3,'padding','circular'); % center about 0
    for iSeason = 1:4
        polarplot((theta(seasonCell{iSeason})),repmat(nanmean(smoothMedSOL(seasonCell{iSeason})),...
            size(seasonCell{iSeason})),'color',[seasonColors(iSeason,:),0.4],'linewidth',15);
        hold on;
    end
    polarplot(theta,smoothMedSOL,'color','k','linewidth',3);
    
    
    title(titleLabels_sol{iSun});
    rticks(0:20:200);
    pax = gca;
    pax.ThetaTick = linspace(0,360,13);
    pax.ThetaTickLabels = useMonths;
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 14;
    grid on;
    
    subplot(rows,cols,3);
    smoothSleepWakeMed = imgaussfilt(sleepWake_med,3,'padding','circular'); % center about 0
    polarplot(theta,smoothSleepWakeMed,'color',sleepWakeColors(iSun,:),'linewidth',3);
    hold on;
    polarplot(theta,repmat(720,[1,366]),'color','k','linewidth',1);
    hold on;
    title(titleLabels_sol{iSun});
    rticks(720);
    rticklabels([]);
    rlim([720-360 720+160]);
    pax = gca;
    pax.ThetaTick = linspace(0,360,13);
    pax.ThetaTickLabels = useMonths;
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 14;
    grid on;
end