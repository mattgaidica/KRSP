cd '/Users/mattgaidica/Documents/MATLAB/KRSP';

loadPath = '/Users/mattgaidica/Dropbox (University of Michigan)/from_box/KRSP Axy Data/Temp';

sqkey = readtable('sqkey.txt');
longevity = readtable(fullfile(pwd,'R','krsp_longevity.csv'));
cone_counts = readtable('R/krsp_cone_counts.csv');
midden_cones = readtable('R/krsp_midden_cones.csv');

exportPath = '/Users/mattgaidica/Documents/MATLAB/KRSP/export';
seasonShiftDays = 56; %  center by light, *seasonDoys(1) = 311
% seasonShiftDays = 56 - 21; %  center by temp
sIds = round(linspace(1,366,5));
seasonDoys = circshift(1:366,seasonShiftDays);
useDoys = {seasonDoys(sIds(1):sIds(2)),seasonDoys(sIds(2):sIds(3)),...
    seasonDoys(sIds(3):sIds(4)),seasonDoys(sIds(4):sIds(5))};
sTitles = {'Winter','Spring','Summer','Autumn'};
seasonLabels = {'Winter','Spring','Summer','Autumn'};
seasonLookup = {};
seasonIndex = [];
for ii = 1:366
    for jj = 1:4
        if ismember(ii,useDoys{jj})
            seasonLookup{ii} = seasonLabels(jj);
            seasonIndex(ii) = jj;
        end
    end
end
seasonLookup = string(seasonLookup);
clc
warning ('off','all');
if do
    weather = readtable('HainesJunction_DailyTemps_Master.csv');
    Tss = makeTss(2014:2020);
    sq_sex = [];
    sq_is_preg = [];
    sq_odba = [];
    sq_nest = [];
    sq_odba_z = [];
    sq_odba_std = [];
    sq_odba_max = [];
    sq_awake = [];
    sq_sqkeyrow = [];
    sq_ids = [];
    sq_ids_un = [];
    sq_doys = [];
    sq_years = [];
    sq_dayLength = [];
    sq_asleep = [];
    iRow = 0;
    recordingId = 0;
    
    alignBySunrise = false; % false = aligns to t=00:00
    
    overlapStats = [];
    overlapMeta = table;
    iNestCount = 0;
    mean_doys = [];
    sq_inNestMin = [];
    sq_asleepMin = [];
    
    trans_all = [];
    trans_secs = [];
    trans_months = [];
    trans_doys = [];
    trans_type = [];
    
    trans_at = [];
    trans_to = [];
    trans_on = [];
    trans_is = [];
    trans_yr = [];
    trans_pg = [];
    
    y_weather = [];
    y_dayLength = [];
    y_asleep = [];
    y_nest = [];
    x_odba = [];
    sqs_odba = cell(1,366);
    sqs_asleep = cell(1,366);
    sqs_nest = cell(1,366);
    sqs_trans = cell(1,366);
    sqs_asleepDay = cell(1,366);
    sqs_asleepNight = cell(1,366);
    sqs_transDay = cell(1,366);
    sqs_transNight = cell(1,366);
    
    trans_dayRatio_isq = [];
    trans_night_isq = [];
    trans_day_isq = [];
    qb_day_sq = [];
    qb_night_sq = [];
    for iSq = 1:size(sqkey,1)
        T = loadTStruct(iSq,sqkey,Tss);
        if isempty(T) || sqkey.is_preg(iSq)
            continue;
        end
        
        Tawake = make_Tawake(T); % transition table

        % /Users/matt/Documents/MATLAB/KRSP/analyze_circCorrSleep.m
        dtdoys = day(T.datetime,'dayofyear');
        dtyears = year(T.datetime);
        [undoys,IA] = unique(dtdoys);
        unyears = dtyears(IA);
        % !! what if I just use sqkey.squirrel_id(iSq) here?
        recordingId = recordingId + 1;
        for iDoy = 1:numel(undoys)
            theseDoys = find(dtdoys == undoys(iDoy));
            if numel(theseDoys) == 1440 % require full day for now
                if alignBySunrise
                    sunrise = Tss.sunrise(Tss.year == unyears(iDoy) & Tss.doy == undoys(iDoy)); 
                    afterSunrise = Tss.sunset(Tss.year == unyears(iDoy) & Tss.doy == undoys(iDoy));
                    closestId = closest(secDay(T.datetime(theseDoys)),secDay(sunrise)); % center on sunrise
                    theseDoys = theseDoys(closestId) - 720:theseDoys(closestId) + 720-1;
                end
                if min(theseDoys) > 1 && max(theseDoys) < numel(T.datetime)
                    iRow = iRow + 1;
                    sq_sqkeyrow(iRow) = iSq;
                    sq_is_preg(iRow) = sqkey.is_preg(iSq);
                    sq_ids(iRow) = recordingId;
                    sq_ids_un(iRow) = sqkey.squirrel_id(iSq);
                    sq_sex(iRow) = strcmp(sqkey.sex{iSq},'M'); % 0 = Female, 1 = Male
                    sq_nest(iRow,:) = strcmp(T.nest(theseDoys),'Nest'); % 0 = out, 1 = in
                    sq_doys(iRow) = undoys(iDoy);
                    sq_odba(iRow,:) = T.odba(theseDoys);
                    sq_odba_z(iRow,:) = T.odba_z(theseDoys);
                    sq_odba_max(iRow,:) = T.odba_max(theseDoys);
                    sq_years(iRow) = year(T.datetime(theseDoys(1)));
                    sq_dayLength(iRow) = Tss.day_length(Tss.year == unyears(iDoy) & Tss.doy == undoys(iDoy));
                    sq_awake(iRow,:) = T.awake(theseDoys);
                    sq_asleep(iRow,:) = T.asleep(theseDoys);
                end
            end
        end
        
        % /Users/matt/Documents/MATLAB/KRSP/analyze_nestSleepOverlap.m
        T.nest_bin = strcmp(T.nest,'Nest');
        thisNestMin = sum(T.nest_bin) ./ size(T,1); % fraction of day
        if thisNestMin*24 > 1 && isValidT(T,true)% must be in nest for at least an hour
            iNestCount = iNestCount + 1;
            sq_inNestMin(iNestCount) = thisNestMin;
            sq_asleepMin(iNestCount) = sum(T.asleep) ./ size(T,1);

            overlapStats(iNestCount,1) = sum(T.nest_bin & T.awake) / size(T,1); % in-awake
            overlapStats(iNestCount,2) = sum(T.nest_bin & ~T.awake) / size(T,1); % in-asleep
            overlapStats(iNestCount,3) = sum(~T.nest_bin & T.awake) / size(T,1); % out-awake
            overlapStats(iNestCount,4) = sum(~T.nest_bin & ~T.awake) / size(T,1); % out-asleep
            recSessUnDoys = unique(day(T.datetime,'dayofyear'));
            mean_doys(iNestCount) = recSessUnDoys(round(numel(recSessUnDoys)/2));
            
            % new, for Ben
            overlapMeta.isq{iNestCount} = iSq;
            overlapMeta.is_preg{iNestCount} = sqkey.is_preg(iSq);
            overlapMeta.squirrelId{iNestCount} = sqkey.squirrel_id(iSq);
            overlapMeta.in_awake{iNestCount} = overlapStats(iNestCount,1);
            overlapMeta.in_asleep{iNestCount} = overlapStats(iNestCount,2);
            overlapMeta.out_awake{iNestCount} = overlapStats(iNestCount,3);
            overlapMeta.out_asleep{iNestCount} = overlapStats(iNestCount,4);
            overlapMeta.sex{iNestCount} = strcmp(sqkey.sex{iSq},'M'); % male = 1
            overlapMeta.is_mast{iNestCount} = ismember(sqkey.year(iSq),[2014,2019]);
            overlapMeta.year{iNestCount} = sqkey.year(iSq);
            overlapMeta.start_dt{iNestCount} = T.datetime(1);
            overlapMeta.end_dt{iNestCount} = T.datetime(end);
            overlapMeta.durationInMinutes{iNestCount} = minutes(T.datetime(end)-T.datetime(1));
            overlapMeta.inNestInMinutes{iNestCount} = sum(T.nest_bin);
            overlapMeta.asleepInMinutes{iNestCount} = sum(~T.awake);
            % add year, start date, end date, start time, end time
            for iSeason = 1:4
                if ismember(mean_doys(iNestCount),seasonDoys(sIds(iSeason):sIds(iSeason+1)))
                    overlapMeta.meanSeason{iNestCount} = iSeason;
                end
            end
        end

        % /Users/matt/Documents/MATLAB/KRSP/analyze_sleepTransitionsAtNight.m
        trans_at = [trans_at secDay(Tawake.datetime)'];
        trans_to = [trans_to Tawake.awake'];
        trans_on = [trans_on day(Tawake.datetime,'dayofyear')'];
        trans_is = [trans_is repmat(iSq,size(Tawake.awake'))];
        trans_yr = [trans_yr year(Tawake.datetime)'];
        trans_pg = [trans_pg repmat(sqkey.is_preg(iSq),size(Tawake.awake'))];
        

        % /Users/matt/Documents/MATLAB/KRSP/results_stats.m
        trans_ratio = [];
        trans_night = [];
        trans_day = [];
        qb_day = [];
        qb_night = [];
        for iDoy = 1:numel(undoys)
            sunrise = mean(secDay(Tss.sunrise(Tss.doy == undoys(iDoy))),1);
            sunset = mean(secDay(Tss.sunset(Tss.doy == undoys(iDoy))),1);
            theseDoys = find(dtdoys == undoys(iDoy));
            dayDoys = theseDoys(secDay(T.datetime(theseDoys)) > sunrise & secDay(T.datetime(theseDoys)) < sunset);
            nightDoys = theseDoys(secDay(T.datetime(theseDoys)) < sunrise | secDay(T.datetime(theseDoys)) > sunset);
            [Y,M,D] = datevec(T.datetime(theseDoys(1)));
            useId = find(weather.Date == datetime(Y,M,D));

             % FOR MAST NMAST TABLES
            [Y,M,D] = datevec(T.datetime(theseDoys(1))); %% ??
            if ismember(Y,2014:2020) % all: 2014:2020, mast: [2014,2019], nmast: [2015:2018,2020]
                if numel(theseDoys) == 1440 && ~isempty(useId) && ~isnan(weather.Mean_Temp(useId))
                    y_weather = [y_weather weather.Mean_Temp(useId)];
                    y_dayLength = [y_dayLength mean(Tss.day_length(Tss.doy == undoys(iDoy)),1)];
                    y_asleep = [y_asleep mean(T.asleep(theseDoys))];
                    y_nest = [y_nest sum(strcmp(T.nest(theseDoys),'Nest'))];
                    x_odba = [x_odba mean(T.odba_z(theseDoys))];
                    sqs_odba{undoys(iDoy)} = [sqs_odba{undoys(iDoy)} mean(T.odba(theseDoys))];
                    sqs_asleep{undoys(iDoy)} = [sqs_asleep{undoys(iDoy)} sum(T.asleep(theseDoys))];
                    sqs_asleepDay{undoys(iDoy)} = [sqs_asleepDay{undoys(iDoy)} sum(T.asleep(dayDoys))];
                    sqs_asleepNight{undoys(iDoy)} = [sqs_asleepNight{undoys(iDoy)} sum(T.asleep(nightDoys))];

                    sqs_nest{undoys(iDoy)} = [sqs_nest{undoys(iDoy)} sum(strcmp(T.nest(theseDoys),'Nest'))];

                    sqs_trans{undoys(iDoy)} = [sqs_trans{undoys(iDoy)} sum(abs(diff(T.asleep(theseDoys))))];
                    sqs_transDay{undoys(iDoy)} = [sqs_transDay{undoys(iDoy)} sum(abs(diff(T.asleep(dayDoys))))];
                    sqs_transNight{undoys(iDoy)} = [sqs_transNight{undoys(iDoy)} sum(abs(diff(T.asleep(nightDoys))))];
                end
            end
            qb_day = [qb_day sum(T.asleep(dayDoys)) / numel(theseDoys)];
            qb_night = [qb_night sum(T.asleep(nightDoys)) / numel(theseDoys)];
            trans_ratio = [trans_ratio sum(abs(diff(T.asleep(dayDoys)))) ./ sum(abs(diff(T.asleep(theseDoys))))];
            trans_night = [trans_night sum(abs(diff(T.asleep(nightDoys))))];
            trans_day = [trans_day sum(abs(diff(T.asleep(dayDoys))))];
        end
        qb_day_sq(iSq) = mean(qb_day);
        qb_night_sq(iSq) = mean(qb_night);
        trans_dayRatio_isq(iSq) = mean(trans_ratio);
        trans_night_isq(iSq) = mean(trans_night);
        trans_day_isq(iSq) = mean(trans_day);
    end
    do = false;
    chime;
    writetable(overlapMeta,fullfile('R','nestAsleepOverlap_v3.csv'));
end
warning ('off','all');

%%
doSave = true;
close all

% find unique animals based on rec duration (overlapMeta.duration (in hours))
clc
allIds = [overlapMeta.squirrelId{:}];
[C,IA,IC] = unique(allIds);
recStatsF = [];
recStatsM = [];
for ii = 1:numel(C)
    theseIds = find(C(ii) == allIds);
    fprintf("%i - ",allIds(theseIds(1)));
    if overlapMeta.is_female{theseIds(1)}
        recStatsF = [recStatsF numel(theseIds)];
    else
        recStatsM = [recStatsM numel(theseIds)];
    end
    for jj = 1:numel(theseIds)
        if jj > 1
            fprintf(" | ");
        end
        fprintf("%i/%i:%ihrs",month(overlapMeta.start_dt{theseIds(jj)}),year(overlapMeta.start_dt{theseIds(jj)}),round(overlapMeta.durationInMinutes{theseIds(jj)}/3600));
    end
    fprintf('\n');
end

subplotMargins = [.1,.15]; % [vert, horz]
colors = lines(7);
h = ff(500,600);
rows = 2;
cols = 2;

% histogram by year, color mast season
subplot_tight(rows,cols,1,subplotMargins);
fcounts = histcounts([overlapMeta.year{[overlapMeta.is_female{:}]}],2013.5:2019.5);
mcounts = histcounts([overlapMeta.year{~[overlapMeta.is_female{:}]}],2013.5:2019.5);
b = bar(2014:2019,[mcounts' fcounts'],'stacked');
b(1).FaceColor = colors(6,:);
b(2).FaceColor = colors(7,:);
set(gca,'fontsize',14);
legend({'Male','Female'});
xtickangle(45);
ylabel('Recording Sessions');
grid on;
legend box off;

% histogram by rec sessions
subplot_tight(rows,cols,2,subplotMargins);
fcounts = histcounts(recStatsF,0.5:8.5);
mcounts = histcounts(recStatsM,0.5:8.5);
b = bar(1:8,[mcounts' fcounts'],'stacked');
b(1).FaceColor = colors(6,:);
b(2).FaceColor = colors(7,:);
set(gca,'fontsize',14);
legend({'Male','Female'});
ylabel('Recording Sessions');
xlabel({'Sessions per Squirrel'});
grid on;
legend box off;

% timeline by recording
subplot_tight(rows,cols,[3,4],subplotMargins);
[~,k] = sort(day([overlapMeta.start_dt{:}],'dayofyear'));

useShift = 1:366;
for ii = 1:size(overlapMeta,1)
    shiftDoy = day(overlapMeta.start_dt{k(ii)},'dayofyear');
    recDays = ceil(overlapMeta.durationInMinutes{ii} / 3600);
    xs = circshift(useShift,-shiftDoy);
    xs = xs(1:recDays);
    if overlapMeta.is_female{ii}
        useColor = colors(7,:);
    else
        useColor = colors(6,:);
    end

    plot(xs,repmat(ii,size(xs)),'.','color',useColor); % add color M/F
    hold on;
end
lns = [];
lns(1) = plot(0,0,'-','linewidth',3,'color',colors(6,:));
lns(2) = plot(0,0,'-','linewidth',3,'color',colors(7,:));

xlim([1,366]);
ylim([1,size(overlapMeta,1)]);
yticks(ylim);
grid on;
set(gca,'fontsize',14);
ylabel('Recording Session');
xlabel('Day of Year');
legend(lns,{'Male','Female'},'location','southeast');
legend box off;

XY = [-.25,1.1;-.25,1.1;-.1,1.1];
addFigureLabels(h,XY);
set(h,'PaperPositionMode','auto');
if doSave
    saveas(h,fullfile(exportPath,'dataLandscape.eps'),'epsc');
    saveas(h,fullfile(exportPath,'dataLandscape.jpg'),'jpg');
    close(h);
end
%% find ideal season doy
close all
sTitles = {'winter','spring','summer','autumn'};
sIds = round(linspace(1,366,5));
all_diff = [];
iNestCount = 0;
useshift = 40:80;
for ii = useshift
    disp(ii);
    seasonDoys = circshift(1:366,ii);
    meanDayLength_spring = mean(Tss.day_length(seasonDoys(sIds(2):sIds(3))));
    meanDayLength_autumn = mean(Tss.day_length(seasonDoys(sIds(4):sIds(5))));
    iNestCount = iNestCount + 1;
    all_diff(iNestCount) = abs(meanDayLength_autumn - meanDayLength_spring);
end
ff(1200,400);
plot(all_diff);
[v,k] = min(all_diff);
title(sprintf('min shift has discrepency of %1.2fs at %1.0f',v,useshift(k)));