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
% close all
useMonths =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
titleLabels_sol = {'Sunrise SOL','Sunset SOL'};
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
        useIds = find(T_SOL.isSunrise == sunriseKey(iSun) & ismember(T_SOL.doy,useDoys) & ismember(T_SOL.is_mast,[0,1]));
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
    title('Sleep/Wake Timing');
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
%% how consistent are SOLs for an individual rec session?
SOLs = [];
curSq = NaN;
unSqs = unique(T_SOL.sq_key_id);
sunriseKey = [1,0]; % goes: sunrise, sunset

nSurr = 1000;
catArr = [];

for iSurr = 1:nSurr+1 % first entry is real
    fprintf("iSurr %i/%i\n",iSurr,nSurr);
    catCount = 0;
    for iSun = 2
        for iSq = 1:numel(unSqs)
            useIds = find(T_SOL.sq_key_id == unSqs(iSq) & T_SOL.isSunrise == sunriseKey(iSun));
            if iSurr > 1 % overwrite with surrogates
                surrIds = find(T_SOL.isSunrise == sunriseKey(iSun));
                useIds = surrIds(randperm(numel(surrIds),numel(useIds)));
            end
            if numel(useIds) > 2
                for ii = 2:numel(useIds)-1
                    catCount = catCount + 1;
                    dayBefore = T_SOL.SOL(useIds(ii-1));
                    curDay = T_SOL.SOL(useIds(ii));
                    dayAfter = T_SOL.SOL(useIds(ii+1));
                    if dayBefore < curDay && curDay > dayAfter
                        catArr(iSurr,catCount) = 1; %#ok<*SAGROW>
                    end
                    if dayBefore > curDay && curDay < dayAfter
                        catArr(iSurr,catCount) = 1;
                    end
                    if dayBefore < curDay && curDay < dayAfter
                        catArr(iSurr,catCount) = 2;
                    end
                    if dayBefore > curDay && curDay > dayAfter
                        catArr(iSurr,catCount) = 2;
                    end
                    if dayBefore > curDay && curDay == dayAfter
                        catArr(iSurr,catCount) = 3;
                    end
                    if dayBefore < curDay && curDay == dayAfter
                        catArr(iSurr,catCount) = 3;
                    end
                    if dayBefore == curDay && curDay < dayAfter
                        catArr(iSurr,catCount) = 3;
                    end
                    if dayBefore == curDay && curDay > dayAfter
                        catArr(iSurr,catCount) = 3;
                    end
                    if dayBefore == curDay && curDay == dayAfter
                        catArr(iSurr,catCount) = 3;
                    end
                end
            end
        end
    end
end
surrArr = catArr(2:end,:);
catArr = catArr(1,:);

%% ^ analyze, result: it actually seems that dayon-dayoff-dayon is less likely in individuals (rest day hypothesis)
nCats = numel(unique(catArr));
pArr = [];
avgArr = [];
bins = 0.5:1:nCats+0.5;
catHistCounts = histcounts(catArr,bins);
for iSurr = 1:nSurr
    surrHistCounts = histcounts(surrArr(iSurr,:),bins);
    avgArr(iSurr,:) = surrHistCounts;
% %     for iCat = 1:nCats
% %         if catHistCounts(iCat) > surrHistCounts(iCat)
% %             pArr(iSurr,iCat) = 1;
% %         elseif catHistCounts(iCat) < surrHistCounts(iCat)
% %             pArr(iSurr,iCat) = -1;
% %         else
% %             pArr(iSurr,iCat) = 0;
% %         end
% %     end
end

close all;
ff(1000,300);
for iCat = 1:nCats
    subplot(1,3,iCat);
    histogram(avgArr(:,iCat));
    xline(catHistCounts(iCat),'linewidth',3,'color','k');
    xline(prctile(avgArr(:,iCat),5),'r');
    xline(prctile(avgArr(:,iCat),95),'r');
    title(sprintf("Cat %i",iCat));
end

%% do individuals have a consistent SOL? result: seems that high SOL at night associate with large SOL std (crappy sleepers?)
meanArr = [];
stdArr = [];
arrCount = 0;

sunriseKey = [1,0]; % goes: sunrise, sunset
colors = lines(3);
for iSun = 1:2
    for iSq = 1:numel(unSqs)
        useIds = find(T_SOL.sq_key_id == unSqs(iSq) & T_SOL.isSunrise == sunriseKey(iSun));
        meanArr(iSq,iSun) = mean(T_SOL.SOL(useIds));
        stdArr(iSq,iSun) = std(T_SOL.SOL(useIds));
    end
end
% [~,k] = sort(meanArr(:,2));
[~,k] = sort(EEArr(:,2));
close all
ff(1200,400);
subplot(121);
plot(meanArr(k,1),'color',colors(3,:),'linewidth',2);
hold on;
plot(meanArr(k,2),'color',repmat(0.2,[3,1]),'linewidth',2);
xlim([1 size(meanArr,1)]);
ylim([0 200]);

subplot(122);
plot(stdArr(k,1),'color',colors(3,:),'linewidth',2);
hold on;
plot(stdArr(k,2),'color',repmat(0.2,[3,1]),'linewidth',2);
xlim([1 size(meanArr,1)]);


%%
T_SOL_sunrise = T_SOL(T_SOL.isSunrise==1,:);
lme = fitlme(T_SOL_sunrise,'SOL~doy+season*is_mast+sex+(1|sq_key_id)')

%% do squirrels actually sleep if they go into the nest early?
seasonColors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
sunriseKey = [1,0]; % goes: sunrise, sunset
EEArr = [];
SOLArr = [];
colorsArr = [];
recCount = 0;
for iSun = 1:2
    useIds = find(T_SOL.isSunrise == sunriseKey(iSun));
    for ii = 1:numel(useIds)
        recCount = recCount + 1;
        EEArr(recCount,iSun) = T_SOL.firstEE(useIds(ii));
        SOLArr(recCount,iSun) = T_SOL.SOL(useIds(ii));
        colorsArr(recCount,:) = seasonColors(T_SOL.season(useIds(ii),:),:);
    end
end
clc
ff(1200,400);
subplot(121);
scatter(EEArr(:,1)-720,SOLArr(:,1),5,colorsArr,'filled');
[r,p] = corr(EEArr(:,1),SOLArr(:,1));
fprintf("Sunrise Exit x Sleep Offset Latency: r = %1.2f, p = %1.2e\n",r,p);
subplot(122);
scatter(EEArr(:,2)-720,SOLArr(:,2),5,colorsArr,'filled');
[r,p] = corr(EEArr(:,2),SOLArr(:,2));
fprintf("Sunset Exit x Sleep Onset Latency: r = %1.2f, p = %1.2e\n",r,p);

%% find outlier entries
seasonColors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
sunriseKey = [1,0]; % goes: sunrise, sunset
EEArr = [];
SOLArr = [];
colorsArr = [];
sqKeyArr = [];
recCount = 0;
useIds = find(T_SOL.isSunrise == 1 & T_SOL.firstEE-720 < -300);
for ii = 1:numel(useIds)
    recCount = recCount + 1;
    EEArr(recCount) = T_SOL.firstEE(useIds(ii));
    SOLArr(recCount) = T_SOL.SOL(useIds(ii));
    colorsArr(recCount,:) = seasonColors(T_SOL.season(useIds(ii),:),:);
    sqKeyArr(recCount) = T_SOL.sq_key_id(useIds(ii));
end
% % close all
% % ff(600,400);
% % scatter(EEArr-720,SOLArr,5,colorsArr,'filled');
% % text(EEArr-720,SOLArr,compose('%i',sqKeyArr));

clc
sqKeyArr_un = unique(sqKeyArr);
fracBad = [];
for ii = 1:numel(sqKeyArr_un)
    useIds = find(T_SOL.isSunrise == 1 & T_SOL.sq_key_id == sqKeyArr_un(ii));
    eeArr = [];
    for jj = 1:numel(useIds)
        eeArr(jj) = T_SOL.firstEE(useIds(jj))-720;
    end
    fracBad(ii) = sum(eeArr < -200) / numel(eeArr);
    fprintf("sqKeyId %i: %i/%i\n",sqKeyArr_un(ii),sum(eeArr < -200),numel(eeArr));
end
badIds = sqKeyArr_un(fracBad > .5);
fprintf("bad: %s\n",compose("%i",badIds));
% bad: 392
% bad: 393
% bad: 394
% bad: 396
% bad: 398
% bad: 401
% bad: 402
% bad: 404
% bad: 406
% bad: 407
% bad: 409
% bad: 411
% bad: 412
% bad: 415
% bad: 416
% bad: 419
% bad: 420
% bad: 422