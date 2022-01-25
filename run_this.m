%% how many records have files? This is really the starting point
iFile = 0;
for iSq = 1:size(sqkey,1)
    if isempty(sqkey.filename{iSq})
        iFile = iFile + 1;
    end
end
fprintf("%i no file\n",iFile);

%%
iValid = 0;
iFemalePreg = 0;
femaleRecDays = [];
femaleSqs = [];
for iSq = 1:size(sqkey,1)
    if isempty(sqkey.filename{iSq}) || sqkey.isValid(iSq) == 0
        continue;
    end
    T = loadTStruct(iSq,sqkey,Tss);
    if isempty(T)
        continue;
    end
    iValid = iValid + 1;

    if strcmp(sqkey.sex_status{iSq},'lactating') || strcmp(sqkey.sex_status{iSq},'pregnant') ||...
            strcmp(sqkey.sex_status{iSq},'Pre-pregnancy')
        iFemalePreg = iFemalePreg + 1;
        femaleRecDays(iFemalePreg) = size(T,1)/1440;
        femaleSqs(iFemalePreg) = iSq;
    end
end

%% ^goes with, for females
uniqFids = unique(sqkey.squirrel_id(femaleSqs));
fyears = sqkey.year(femaleSqs);

fprintf("%i valid\n",iValid);
fprintf("%i rec (%i uniq) female preg\n",iFemalePreg,numel(uniqFids));

allConds = lower(sqkey.treatment(femaleSqs));
unConds = unique(allConds);
histArr = [];
for ii = 1:numel(unConds)
    histArr(ii) = sum(strcmp(allConds,unConds{ii}));
end

colors = lines(1);
close all;
ff(800,300);
subplot(131);
histogram(femaleRecDays,10,'facecolor',colors(1,:),'facealpha',1);
xlabel('rec. days in session');
ylabel('rec. sessions');
ylim([0 30]);
title(sprintf('%i rec. sessions\n(%i unique squirrels)',iFemalePreg,numel(uniqFids)));

subplot(132);
histogram(fyears,'facecolor',colors(1,:),'facealpha',1);
ylabel('rec. sessions');
xtickangle(45);
xticks(2014:2020);
ylim([0 30]);

subplot(133);
bar(histArr);
xticklabels(unConds);
xtickangle(45);
ylabel('rec. sessions');
ylim([0 30]);

saveas(gcf,fullfile(exportPath,'female_recStats.png'));
%% preg outline after RITable is made in xcorrAsleep.m
seasonStr = {'Winter','Spring','Summer','Autumn'};
mastStr = {'nMast','Mast'};
pregStr = {'nPreg','Preg'};
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
labelStr = {};
iCount = 0;
barMeanRI = [];
barMeanQB = [];
barStdRI = [];
barStdQB = [];
barColors = [];
isPreg = [];
for iSeason = 1:4
    for iMast = 0:1
        for iPreg = 0:1
            ids = find(RITable.season == iSeason & RITable.is_mast == iMast & RITable.is_preg == iPreg);
            if ~isempty(ids)
                iCount = iCount + 1;
                labelStr{iCount} = sprintf('%s-%s-%s',seasonStr{iSeason},mastStr{iMast+1},pregStr{iPreg+1});
                barMeanRI(iCount) = mean(RITable.RI(ids));
                barStdRI(iCount) = std(RITable.RI(ids));
                barMeanQB(iCount) = mean(RITable.qb(ids));
                barStdQB(iCount) = std(RITable.qb(ids));
                if iMast == 0
                    barColors(iCount,:) = colors(iSeason,:);
                else
                    barColors(iCount,:) = colors(iSeason,:).^2;
                end
                isPreg(iCount) = iPreg;
            end
        end
    end
end

close all;
ff(900,800);
useMean = {barMeanRI,barMeanQB};
useStd = {barStdRI,barStdQB};
axLabels = {'RI (raw value)','QB (Z-score)'};
useOffset = [0.05,0.25];
for iPlot = 1:2
    subplot(2,1,iPlot);
    for ii = 1:size(barColors,1)
        errorbar(ii,useMean{iPlot}(ii),useStd{iPlot}(ii),'linewidth',4,'color',barColors(ii,:));
        hold on;
        plot(ii,useMean{iPlot}(ii),'.','color',barColors(ii,:),'markersize',30);
        if isPreg(ii) == 1
            text(ii,useMean{iPlot}(ii)+useStd{iPlot}(ii)+useOffset(iPlot),'P',...
                'fontsize',16,'color','k','horizontalalignment','center');
        end
    end
    xticklabels(labelStr)
    xticks(1:size(barColors,1));
    ylabel(axLabels{iPlot});
    xtickangle(30);
    xlim([0,size(barColors,1)+1]);
    colororder(barColors);
    title('Mast years are darker, pregnant data marked P');
end
%%
% litter = readtable('litter.csv');
% RITable = readtable(fullfile('R','RITable.csv'));
close all
ff(750,300);
lw1 = 1;
lw2 = 5;
ms1 = 50;
ms2 = 15;
y1 = 0.65;
y2 = 0.70;
lns = [];

rowIds = find(ismember(year(litter.fieldBDate),[2015:2018,2020]));
allnMastDoys = day(litter.fieldBDate(rowIds),'doy');

nmastIds = find(RITable.is_mast == 0);
lns(1) = plot(RITable.doy(nmastIds),RITable.RI(nmastIds),'k.','markersize',ms2);
hold on;
x1 = min(allnMastDoys);
x2 = max(allnMastDoys);
lns(3) = plot([x1,x2],[y1,y1],'k-','linewidth',lw1);
x1 = mean(allnMastDoys) - std(allnMastDoys);
x2 = mean(allnMastDoys) + std(allnMastDoys);
lns(4) = plot([x1,x2],[y1,y1],'k-','linewidth',lw2);
lns(5) = plot(mean(allnMastDoys),y1,'k.','markersize',ms1);
lns(6) = plot(mean(allnMastDoys)-38,y1,'kx','markersize',ms2);

rowIds = find(ismember(year(litter.fieldBDate),[2014,2019]));
allMastDoys = day(litter.fieldBDate(rowIds),'doy');

mastIds = find(RITable.is_mast == 1);
lns(2) = plot(RITable.doy(mastIds),RITable.RI(mastIds),'r.','markersize',ms2);
x1 = min(allMastDoys);
x2 = max(allMastDoys);
plot([x1,x2],[y2,y2],'r-','linewidth',lw1);
x1 = mean(allMastDoys) - std(allMastDoys);
x2 = mean(allMastDoys) + std(allMastDoys);
plot([x1,x2],[y2,y2],'r-','linewidth',lw2);
plot(mean(allMastDoys),y2,'r.','markersize',ms1);
plot(mean(allMastDoys)-38,y2,'rx','markersize',ms2);

legend(lns,{'Non-mast','Mast','Min-max Birth','1 Std. Birth','Mean Birth','Mean-38d'});
xticks([35,127,219,311]);
xticklabels({'Spring→','Summer→','Autumn→','Winter→'});
xlabel('DOY');
ylabel('RI');
title('RI vs. Day of Year');
xlim([1 366]);
set(gca,'fontsize',14);
ylim([0 0.75]);
grid on

%%
close all
ff(600,600);
for iSeason = 1:4
    subplot(2,2,iSeason)
    ids = find(RITable.season == iSeason & RITable.is_mast == 0 & ~isnan(RITable.traps_rec));
    plot(RITable.traps_rec(ids),RITable.qb(ids),'kx');
    if ~isempty(ids)
        [r,p] = corr(RITable.traps_rec(ids),RITable.qb(ids));
        fprintf('%i-nm: r=%1.4f, p=%1.4f\n',iSeason,r,p);
    end
    hold on;
    ids = find(RITable.season == iSeason & RITable.is_mast == 1 & ~isnan(RITable.traps_rec));
    plot(RITable.traps_rec(ids),RITable.qb(ids),'ro');
    if ~isempty(ids)
        [r,p] = corr(RITable.traps_rec(ids),RITable.qb(ids));
        fprintf('%i-ma: r=%1.4f, p=%1.4f\n',iSeason,r,p);
    end
    title(sprintf('Season %i',iSeason));
    xlim([-2.5 2.5]);
    ylim(xlim);
end
%%
all_ids = find(ismember(RITable.season,4) & ~isnan(RITable.traps_rec));
% nmast_ids = find(ismember(RITable.season,3) & RITable.is_mast == 0 & ~isnan(RITable.traps_rec));
useLims = [-2.5 2.5];

close all
ff(900,300);
subplot(131);
x = RITable.traps_rec(all_ids);
y = normalize(RITable.qb(all_ids));
plot(x,y,'k.','markersize',40);
[r,p] = corr(x,y);

f = fit(x,y,'poly1');
hold on;
plot(useLims,f(useLims),'r-');

xticks([useLims(1),0,useLims(2)]);
yticks([useLims(1),0,useLims(2)]);
xlim(useLims);
ylim(useLims);
set(gca,'fontsize',14);
xlabel('Trapping Incidence (Z, TI)');
ylabel('Quiescent Behavior (Z, QB)');
title(sprintf("Summer QB vs. TI\nr = %1.2f, p = %1.2e",r,p));
grid on;

subplot(132);
y = normalize(RITable.qb_day(all_ids));
plot(x,y,'k.','markersize',40);
[r,p] = corr(x,y);

f = fit(x,y,'poly1');
hold on;
plot(useLims,f(useLims),'r-');

xticks([useLims(1),0,useLims(2)]);
yticks([useLims(1),0,useLims(2)]);
xlim(useLims);
ylim(useLims);
set(gca,'fontsize',14);
xlabel('Trapping Incidence (Z, TI)');
ylabel('Quiescent Behavior (Z, QB)');
title(sprintf("Summer Day-QB vs. TI\nr = %1.2f, p = %1.2e",r,p));
grid on;

subplot(133);
y = normalize(RITable.qb_night(all_ids));
plot(x,y,'k.','markersize',40);
[r,p] = corr(x,y);

f = fit(x,y,'poly1');
hold on;
plot(useLims,f(useLims),'r-');

xticks([useLims(1),0,useLims(2)]);
yticks([useLims(1),0,useLims(2)]);
xlim(useLims);
ylim(useLims);
set(gca,'fontsize',14);
xlabel('Trapping Incidence (Z, TI)');
ylabel('Quiescent Behavior (Z, QB)');
title(sprintf("Summer Night-QB vs. TI\nr = %1.2f, p = %1.2e",r,p));
grid on;