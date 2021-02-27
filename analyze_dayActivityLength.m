% setup with /Users/matt/Documents/MATLAB/KRSP/Figures/figure_rhythmicSleep_noDiff.m
% these analyses compare ODBA with temperature, day length, and
% also contains cursory analyses on sleep transitions, used in
% my weekly update
doAnalysis = true;

weatherPath = '/Users/matt/Documents/Data/KRSP/HainesJunction_DailyTemps_Master.csv';
T_weather = readtable(weatherPath);
ssPath = '/Users/matt/Documents/Data/KRSP/SunriseSunset';
files = dir(fullfile(ssPath,'*.txt'));
T_ss = readtable(fullfile(ssPath,files(4).name)); % 366 day year (simplify for now)

%% What is the relationship between day length and activity?
% hypothesis: no correlation, activity per unit sunlight is equal across all doys

cmap_season = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png',366);
dayOdba = [];
nightOdba = [];
colors_season = [];
cLines = lines(7);
colors_lines = [cLines(7,:);cLines(1,:);cLines(3,:);cLines(4,:)];
colors_sex = [];
colors_mast = [];
temps = [];
asleep = [];
asleep_trans = [];
asleep_trans_day = [];

for iSq = 1:numel(sq_ids)
    thisDoy = sq_doys(iSq);
    dayOdba(iSq) = sum(sq_odba_std(iSq,minDay(T_ss.sunrise(thisDoy)):minDay(T_ss.sunset(thisDoy))));
    nightOdba(iSq) = sum([sq_odba_std(iSq,1:minDay(T_ss.sunrise(thisDoy))),...
        sq_odba_z(iSq,minDay(T_ss.sunset(thisDoy)):end)]);
    
    colors_season(iSq,:) = cmap_season(thisDoy,:);
    colors_sex(iSq,:) = colors_lines(sq_sex(iSq)+1,:);
    temps(iSq) = T_weather.Mean_Temp(T_weather.Year == 2016 & T_weather.Julian_Date == thisDoy);
    asleep(iSq) = sum(sq_asleep(iSq,:));
    asleep_trans(iSq) = sum(abs(diff(sq_asleep(iSq,:))) == 1) + rand;
    asleep_trans_night(iSq) = ((sum(abs(diff(sq_asleep(iSq,1:minDay(T_ss.sunrise(thisDoy))))) == 1)...
        + sum(abs(diff(sq_asleep(iSq,minDay(T_ss.sunset(thisDoy)):end))) == 1)) + rand) ...
        ./ (86400 - sq_dayLength(iSq));
    asleep_trans_day(iSq) = (sum(abs(diff(sq_asleep(iSq,minDay(T_ss.sunrise(thisDoy)):minDay(T_ss.sunset(thisDoy))))) == 1)...
        + rand) ./ (sq_dayLength(iSq));
    colors_mast(iSq,:) = colors_lines(ismember(sq_years(iSq),[2014,2019])+3,:); % 2014,2019
end
asleep_trans_night = normalize(asleep_trans_night,'range');
asleep_trans_day = normalize(asleep_trans_day,'range');

%% FIG 1
close all

rows = 4;
cols = 2;
S = 5;
ff(750,1000);

% sleep variance

subplot(rows,cols,1);
scatter(sq_dayLength/60/60,dayOdba,S,colors_season,'filled');
xlabel('dayLength (hrs)');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,2);
scatter(sq_dayLength/60/60,dayOdba,S,colors_sex,'filled');
xlabel('dayLength (hrs)');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
hold on;
mf_lns = []; % do once
mf_lns(1) = plot(nan,nan,'.','markersize',25,'color',colors_lines(1,:));
mf_lns(2) = plot(nan,nan,'.','markersize',25,'color',colors_lines(2,:));
legend(mf_lns,{'Female','Male'},'location','northwest');
title('');

subplot(rows,cols,3);
scatter(temps,dayOdba,S,colors_season,'filled');
xlabel('temp (C)');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,4);
scatter(temps,dayOdba,S,colors_sex,'filled');
xlabel('temp (C)');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,5);
scatter(sq_dayLength/60/60,asleep,S,colors_season,'filled');
xlabel('dayLength (hrs)');
ylabel('asleep (min)');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,6);
scatter(sq_dayLength/60/60,asleep,S,colors_sex,'filled');
xlabel('dayLength (hrs)');
ylabel('asleep (min)');
set(gca,'fontsize',14);
hold on;
title('');

subplot(rows,cols,[rows*cols-1 rows*cols]);
imshow(imresize(imread('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png'),[24*3,366],'nearest'));
hold on;
plot(1:366,(T_ss.day_length./60/60)*3,'w-','linewidth',2);
axis on;
xticks(linspace(min(xlim),max(xlim),13));
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
yticks([0:6:23]*3);
yticklabels(compose('%1.2f',yticks/3));
set(gca,'fontsize',14);
set(gca,'ydir','normal');
grid on;
ylabel('daylight (hrs)');

%% FIG 2
rows = 5;
cols = 2;
S = 5;
ff(750,1000);

subplot(rows,cols,1);
scatter(sq_dayLength/60/60,dayOdba,S,colors_season,'filled');
xlabel('dayLength (hrs)');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,2);
scatter(sq_dayLength/60/60,dayOdba,S,colors_mast,'filled');
xlabel('dayLength (hrs)');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
hold on;
mf_lns = []; % do once
mf_lns(1) = plot(nan,nan,'.','markersize',25,'color',colors_lines(3,:));
mf_lns(2) = plot(nan,nan,'.','markersize',25,'color',colors_lines(4,:));
legend(mf_lns,{'Non-mast','Mast'},'location','northwest');
title('');

subplot(rows,cols,3);
scatter(sq_dayLength/60/60,nightOdba,S,colors_season,'filled');
xlabel('dayLength (hrs)');
ylabel('sum night \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,4);
scatter(sq_dayLength/60/60,nightOdba,S,colors_mast,'filled');
xlabel('dayLength (hrs)');
ylabel('sum night \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,5);
scatter(temps,dayOdba,S,colors_season,'filled');
xlabel('temp (C)');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,6);
scatter(temps,dayOdba,S,colors_mast,'filled');
xlabel('temp (C)');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,7);
scatter(sq_dayLength/60/60,asleep,S,colors_season,'filled');
xlabel('dayLength (hrs)');
ylabel('asleep (min)');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,8);
scatter(sq_dayLength/60/60,asleep,S,colors_mast,'filled');
xlabel('dayLength (hrs)');
ylabel('asleep (min)');
set(gca,'fontsize',14);
hold on;
title('');

subplot(rows,cols,[rows*cols-1 rows*cols]);
imshow(imresize(imread('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png'),[24*3,366],'nearest'));
hold on;
plot(1:366,(T_ss.day_length./60/60)*3,'w-','linewidth',2);
axis on;
xticks(linspace(min(xlim),max(xlim),13));
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
yticks([0:6:23]*3);
yticklabels(compose('%1.2f',yticks/3));
set(gca,'fontsize',14);
set(gca,'ydir','normal');
grid on;
ylabel('daylight (hrs)');

%% FIG 3
rows = 5;
cols = 2;
S = 7;
ff(750,1000);

subplot(rows,cols,1);
scatter(sq_dayLength/60/60,asleep,S,colors_season,'filled');
xlabel('dayLength (hrs)');
ylabel('asleep (min)');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,2);
scatter(sq_dayLength/60/60,asleep,S,colors_sex,'filled');
xlabel('dayLength (hrs)');
ylabel('asleep (min)');
set(gca,'fontsize',14);
hold on;
mf_lns = []; % do once
mf_lns(1) = plot(nan,nan,'.','markersize',25,'color',colors_lines(1,:));
mf_lns(2) = plot(nan,nan,'.','markersize',25,'color',colors_lines(2,:));
legend(mf_lns,{'Female','Male'},'location','northeast');
title('');

subplot(rows,cols,3);
scatter(asleep_trans,dayOdba,S,colors_season,'filled');
xlabel('sleep transitions (all day)');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,4);
scatter(asleep_trans,dayOdba,S,colors_sex,'filled');
xlabel('sleep transitions (all day)');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,5);
scatter(asleep_trans_day,dayOdba,S,colors_season,'filled');
xlabel('sleep transitions (day frac_{norm})');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,6);
scatter(asleep_trans_day,dayOdba,S,colors_sex,'filled');
xlabel('sleep transitions (day frac_{norm})');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,7);
scatter(asleep_trans_night,dayOdba,S,colors_season,'filled');
xlabel('sleep transitions (night frac_{norm})');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,8);
scatter(asleep_trans_night,dayOdba,S,colors_sex,'filled');
xlabel('sleep transitions (night frac_{norm})');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

% Night is similar to day, kind of redundant info
% % subplot(2,2,2);
% % scatter(sq_dayLength/60/60,nightOdba,10,colors_season,'filled');
% % xlabel('dayLength (hrs)');
% % ylabel('sum \DeltaOA');
% % set(gca,'fontsize',14);
% % title('Night');

subplot(rows,cols,[rows*cols-1 rows*cols]);
imshow(imresize(imread('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png'),[24*3,366],'nearest'));
hold on;
plot(1:366,(T_ss.day_length./60/60)*3,'w-','linewidth',2);
axis on;
xticks(linspace(min(xlim),max(xlim),13));
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
yticks([0:6:23]*3);
yticklabels(compose('%1.2f',yticks/3));
set(gca,'fontsize',14);
set(gca,'ydir','normal');
grid on;
ylabel('daylight (hrs)');

%% FIG 4
rows = 5;
cols = 2;
S = 7;
ff(750,1000);

subplot(rows,cols,1);
scatter(sq_dayLength/60/60,asleep,S,colors_season,'filled');
xlabel('dayLength (hrs)');
ylabel('asleep (min)');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,2);
scatter(sq_dayLength/60/60,asleep,S,colors_mast,'filled');
xlabel('dayLength (hrs)');
ylabel('asleep (min)');
set(gca,'fontsize',14);
hold on;
mf_lns = []; % do once
mf_lns(1) = plot(nan,nan,'.','markersize',25,'color',colors_lines(3,:));
mf_lns(2) = plot(nan,nan,'.','markersize',25,'color',colors_lines(4,:));
legend(mf_lns,{'Non-mast','Mast'},'location','northeast');
title('');

subplot(rows,cols,3);
scatter(asleep_trans,dayOdba,S,colors_season,'filled');
xlabel('sleep transitions (all day)');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,4);
scatter(asleep_trans,dayOdba,S,colors_mast,'filled');
xlabel('sleep transitions (all day)');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
hold on;
mf_lns = []; % do once
mf_lns(1) = plot(nan,nan,'.','markersize',25,'color',colors_lines(3,:));
mf_lns(2) = plot(nan,nan,'.','markersize',25,'color',colors_lines(4,:));
legend(mf_lns,{'Non-mast','Mast'},'location','northeast');
title('');

subplot(rows,cols,5);
scatter(asleep_trans_day,dayOdba,S,colors_season,'filled');
xlabel('sleep transitions (day frac_{norm})');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,6);
scatter(asleep_trans_day,dayOdba,S,colors_mast,'filled');
xlabel('sleep transitions (day frac_{norm})');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,7);
scatter(asleep_trans_night,dayOdba,S,colors_season,'filled');
xlabel('sleep transitions (night frac_{norm})');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

subplot(rows,cols,8);
scatter(asleep_trans_night,dayOdba,S,colors_mast,'filled');
xlabel('sleep transitions (night frac_{norm})');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',14);
title('');

% Night is similar to day, kind of redundant info
% % subplot(2,2,2);
% % scatter(sq_dayLength/60/60,nightOdba,10,colors_season,'filled');
% % xlabel('dayLength (hrs)');
% % ylabel('sum \DeltaOA');
% % set(gca,'fontsize',14);
% % title('Night');

subplot(rows,cols,[rows*cols-1 rows*cols]);
imshow(imresize(imread('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png'),[24*3,366],'nearest'));
hold on;
plot(1:366,(T_ss.day_length./60/60)*3,'w-','linewidth',2);
axis on;
xticks(linspace(min(xlim),max(xlim),13));
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
yticks([0:6:23]*3);
yticklabels(compose('%1.2f',yticks/3));
set(gca,'fontsize',14);
set(gca,'ydir','normal');
grid on;
ylabel('daylight (hrs)');

%%
ff(800,900);
subplot(2,cols,1);
scatter(asleep_trans_day,dayOdba,10,colors_season,'filled');
xlabel('sleep transitions (day frac_{norm})');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',22);
title('');

subplot(2,cols,2);
scatter(asleep_trans_night,dayOdba,10,colors_season,'filled');
xlabel('sleep transitions (night frac_{norm})');
ylabel('sum day \DeltaOA');
set(gca,'fontsize',22);
title('');

subplot(2,cols,[2*cols-1 2*cols]);
imshow(imresize(imread('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png'),[24*3,366],'nearest'));
hold on;
plot(1:366,(T_ss.day_length./60/60)*3,'w-','linewidth',2);
axis on;
xticks(linspace(min(xlim),max(xlim),13));
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
yticks([0:6:23]*3);
yticklabels(compose('%1.2f',yticks/3));
set(gca,'fontsize',14);
set(gca,'ydir','normal');
grid on;
ylabel('daylight (hrs)');