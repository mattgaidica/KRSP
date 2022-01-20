%% how many records have files? This is really the starting point
iFile = 0;
for iSq = 1:size(sqkey,1)
    if isempty(sqkey.filename{iSq})
        iFile = iFile + 1;
    end
end
fprintf("%i no file\n",iFile);

%%
longevity = readtable('krsp_longevity.csv');
cone_counts = readtable('krsp_cone_counts.csv');
midden_cones = readtable('krsp_midden_cones.csv');
iValid = 0;
iLon = 0;
iMidden = 0;
iCones = 0;
for iSq = 1:size(sqkey,1)
    if isempty(sqkey.filename{iSq})
        continue;
    end
%     T = loadTStruct(iSq,sqkey,Tss);
    if sqkey.isValid(iSq)
        iValid = iValid + 1;
        lonId = find(longevity.squirrel_id == sqkey.squirrel_id(iSq));
        if ~isempty(lonId)
            iLon = iLon + 1;
        end
        midId = find(midden_cones.squirrel_id == sqkey.squirrel_id(iSq) &...
            midden_cones.year == sqkey.year(iSq));
        if ~isempty(midId)
            iMidden = iMidden + 1;
        end
    end
end
fprintf("%i valid\n",iValid);
fprintf("%i longevity\n",iLon);
fprintf("%i middens\n",iMidden);

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