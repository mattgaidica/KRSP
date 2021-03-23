loadPath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
sqkey = readtable('sqkey.txt');
exportPath = '/Users/matt/Documents/MATLAB/KRSP/export';
clc
if do
    weather = readtable('HainesJunction_DailyTemps_Master.csv');
    Tss = makeTss(2014:2020);
    sq_sex = [];
    sq_odba = [];
    sq_odba_z = [];
    sq_odba_std = [];
    sq_odba_max = [];
    sq_awake = [];
    sq_ids = [];
    sq_doys = [];
    sq_years = [];
    sq_dayLength = [];
    sq_asleep = [];
    iRow = 0;
    squirrelId = 0;
    
    alignBySunrise = false; % false = aligns to t=00:00
    
    overlapStats = [];
    iCount = 0;
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
    for iSq = 1:size(sqkey,1)
        T = loadTStruct(iSq,sqkey,Tss);
        if isempty(T)
            continue;
        end
        
        if strcmp(sqkey.sex_status,'lactating') | strcmp(sqkey.sex_status,'pregnant') | strcmp(sqkey.sex_status,'Pre-pregnancy')
            continue;
        end
        Tawake = make_Tawake(T); % transition table
        
        % /Users/matt/Documents/MATLAB/KRSP/xcorrAsleep.m
%         if numel(T.asleep) >= sq_xcorr_l
%             [c,sq_xcorr_lags] = xcorr(T.asleep,sq_xcorr_l,'coeff');
%             sq_xcorr(size(sq_xcorr,1)+1,:) = c';
%             sq_xcorr_doys = [sq_xcorr_doys;day(T.datetime(round(size(T,1)/2)),'dayofyear')];
%         end

        % /Users/matt/Documents/MATLAB/KRSP/analyze_circCorrSleep.m
        dtdoys = day(T.datetime,'dayofyear');
        dtyears = year(T.datetime);
        [undoys,IA] = unique(dtdoys);
        unyears = dtyears(IA);
        squirrelId = squirrelId + 1;
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
                    sq_ids(iRow) = squirrelId;
                    sq_sex(iRow) = strcmp(sqkey.sex{iSq},'M'); % 0 = Female, 1 = Male
                    sq_doys(iRow) = undoys(iDoy);
                    sq_odba(iRow,:) = T.odba(theseDoys);
                    sq_odba_z(iRow,:) = T.odba_z(theseDoys);
                    %                             sq_odba_std(iRow,:) = T.odba_std(theseDoys);
                    sq_odba_max(iRow,:) = T.odba_max(theseDoys);
                    sq_years(iRow) = year(T.datetime(theseDoys(1)));
                    sq_dayLength(iRow) = Tss.day_length(Tss.year == unyears(iDoy) & Tss.doy == undoys(iDoy));
                    sq_awake(iRow,:) = T.awake(theseDoys);
                    sq_asleep(iRow,:) = T.asleep(theseDoys);
                end
            end
        end
        
        % /Users/matt/Documents/MATLAB/KRSP/analyze_nestSleepOverlap.m
        iCount = iCount + 1;
        T.nest_bin = strcmp(T.nest,'Nest');
        overlapStats(iCount,1) = sum(T.nest_bin & T.awake) / size(T,1); % in-awake
        overlapStats(iCount,2) = sum(T.nest_bin & ~T.awake) / size(T,1); % in-asleep
        overlapStats(iCount,3) = sum(~T.nest_bin & T.awake) / size(T,1); % out-awake
        overlapStats(iCount,4) = sum(~T.nest_bin & ~T.awake) / size(T,1); % out-asleep
        mean_doys(iCount) = mean(unique(day(T.datetime,'dayofyear')));
        sq_inNestMin(iCount) = sum(T.nest_bin) ./ size(T,1);
        sq_asleepMin(iCount) = sum(T.asleep) ./ size(T,1);

        % /Users/matt/Documents/MATLAB/KRSP/analyze_sleepTransitionsAtNight.m
        trans_at = [trans_at secDay(Tawake.datetime)'];
        trans_to = [trans_to Tawake.awake'];
        trans_on = [trans_on day(Tawake.datetime,'dayofyear')'];
        trans_is = [trans_is repmat(iSq,size(Tawake.awake'))];
        trans_yr = [trans_yr year(Tawake.datetime)'];
        
        % /Users/matt/Documents/MATLAB/KRSP/results_stats.m
        for iDoy = 1:numel(undoys)
            sunrise = mean(secDay(Tss.sunrise(Tss.doy == undoys(iDoy))),1);
            sunset = mean(secDay(Tss.sunset(Tss.doy == undoys(iDoy))),1);
            theseDoys = find(dtdoys == undoys(iDoy));
            dayDoys = theseDoys(secDay(T.datetime(theseDoys)) > sunrise & secDay(T.datetime(theseDoys)) < sunset);
            nightDoys = theseDoys(secDay(T.datetime(theseDoys)) < sunrise | secDay(T.datetime(theseDoys)) > sunset);
            [Y,M,D] = datevec(T.datetime(theseDoys(1)));
            useId = find(weather.Date == datetime(Y,M,D));
% %             if ~ismember(Y,[2014,2019])
% %                 continue;
% %             end
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
    end
    do = false;
    chime;
end