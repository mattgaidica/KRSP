weather = readtable('HainesJunction_DailyTemps_Master.csv');
sqkey = readtable('sqkey');
filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
Tss_doys = day(Tss.sunrise,'dayofyear');

unyears = unique(sqkey.year);
unyears = unyears(~isnan(unyears));
sqs_temp = NaN(numel(unyears),366);
for iYear = 1:numel(unyears)
    for iDoy = 1:366
        useId = find(weather.Year == unyears(iYear) & weather.Julian_Date == iDoy);
        if ~isempty(useId)
            sqs_temp(iYear,iDoy) = weather.Mean_Temp(useId);
        end
    end
end

nstd = nanstd(sqs_temp);
v = mean(nanmean(sqs_temp));
fprintf('Mean: %1.2f + %1.2f\n',v,std(nstd));

[v,k] = max(nanmean(sqs_temp));
fprintf('Hottest day: %s, doy=%i, %1.2f + %1.2f\n',datetime(2016,1,1)+k-1,k,v,nstd(k));

[v,k] = min(nanmean(sqs_temp));
fprintf('Coldest day: %s, doy=%i, %1.2f + %1.2f\n',datetime(2016,1,1)+k-1,k,v,nstd(k));

[~,kax] = max(Tss.day_length);
[~,kin] = min(Tss.day_length);
fprintf('maxLength: %1.2f hrs %s doy=%i, %1.2f hrs %s doy=%i\n',...
    max(Tss.day_length)/3600,datetime(2016,1,1)+kax-1,kax,min(Tss.day_length)/3600,datetime(2016,1,1)+kin-1,kin);
fprintf('avg day: %1.2f + %1.2f\n',mean(Tss.day_length)/3600,std(Tss.day_length)/3600);


%% corrs
if do
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
        if ~isempty(sqkey.filename{iSq}) %&& ~any(ismember(sqkey.year(iSq),[2014,2019])) % ~(strcmp(sqkey.source{iSq},'ES') &&
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            T = detect_sleepWake(T);
            dtdoys = day(T.datetime,'dayofyear');
            undoys = unique(dtdoys);
            squirrelId = squirrelId + 1;
            for iDoy = 1:numel(undoys)
                sunrise = secDay(Tss.sunrise(Tss_doys == undoys(iDoy)));
                sunset = secDay(Tss.sunset(Tss_doys == undoys(iDoy)));
                theseDoys = find(dtdoys == undoys(iDoy));
                dayDoys = theseDoys(secDay(T.datetime(theseDoys)) > sunrise & secDay(T.datetime(theseDoys)) < sunset);
                nightDoys = theseDoys(secDay(T.datetime(theseDoys)) < sunrise | secDay(T.datetime(theseDoys)) > sunset);
                [Y,M,D] = datevec(T.datetime(theseDoys(1)));
                useId = find(weather.Date == datetime(Y,M,D));
                if numel(theseDoys) == 1440 && ~isempty(useId) && ~isnan(weather.Mean_Temp(useId))
                    y_weather = [y_weather weather.Mean_Temp(useId)];
                    y_dayLength = [y_dayLength Tss.day_length(undoys(iDoy))];
                    y_asleep = [y_asleep mean(T.asleep(theseDoys))];
                    y_nest = [y_nest sum(strcmp(T.nest(theseDoys),'Nest'))];
                    x_odba = [x_odba mean(T.odba_z(theseDoys))];
                    sqs_odba{undoys(iDoy)} = [sqs_odba{undoys(iDoy)} mean(T.odba_z(theseDoys))];
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
    end
    do = false;
end

%%
figure;
plot(normalize(nanmean(sqs_temp),'scale'));
hold on;
plot(normalize(Tss.day_length,'scale'));
plot(normalize(cellfun(@mean,sqs_odba),'scale'));

%%
clc
sTitles = {'winter','spring','summer','fall'};
sIds = round(linspace(1,366,5));
seasonDoys = circshift(1:366,60);
for iS = 1:4
    meanAsleep = nanmean([sqs_asleep{seasonDoys(sIds(iS):sIds(iS+1))}]);
    stdAsleep = nanstd([sqs_asleep{seasonDoys(sIds(iS):sIds(iS+1))}]);
    fprintf('%s: asleep: %1.2f + %1.2f\n',sTitles{iS},meanAsleep/60,stdAsleep/60);
    
    meanAsleep = nanmean([sqs_asleepDay{seasonDoys(sIds(iS):sIds(iS+1))}]);
    stdAsleep = nanstd([sqs_asleepDay{seasonDoys(sIds(iS):sIds(iS+1))}]);
    fprintf('%s: DAYasleep: %1.2f + %1.2f\n',sTitles{iS},meanAsleep/60,stdAsleep/60);
    
    meanAsleep = nanmean([sqs_asleepNight{seasonDoys(sIds(iS):sIds(iS+1))}]);
    stdAsleep = nanstd([sqs_asleepNight{seasonDoys(sIds(iS):sIds(iS+1))}]);
    fprintf('%s: NIGHTasleep: %1.2f + %1.2f\n\n',sTitles{iS},meanAsleep/60,stdAsleep/60);
    
    meanTrans = nanmean([sqs_trans{seasonDoys(sIds(iS):sIds(iS+1))}]);
    stdTrans = nanstd([sqs_trans{seasonDoys(sIds(iS):sIds(iS+1))}]);
    fprintf('%s: trans: %1.0f + %1.0f\n',sTitles{iS},meanTrans,stdTrans);
    
    meanTrans = nanmean([sqs_transDay{seasonDoys(sIds(iS):sIds(iS+1))}]);
    stdTrans = nanstd([sqs_transDay{seasonDoys(sIds(iS):sIds(iS+1))}]);
    fprintf('%s: DAYtrans: %1.0f + %1.0f\n',sTitles{iS},meanTrans,stdTrans);
    
    meanTrans = nanmean([sqs_transNight{seasonDoys(sIds(iS):sIds(iS+1))}]);
    stdTrans = nanstd([sqs_transNight{seasonDoys(sIds(iS):sIds(iS+1))}]);
    fprintf('%s: NIGHTtrans: %1.0f + %1.0f\n\n',sTitles{iS},meanTrans,stdTrans);
end

meanAsleep = nanmean([sqs_asleep{seasonDoys(1:366)}]);
    stdAsleep = nanstd([sqs_asleep{seasonDoys(1:366)}]);
    fprintf('%s: asleep: %1.2f + %1.2f\n','All',meanAsleep/60,stdAsleep/60);
    
    meanAsleep = nanmean([sqs_asleepDay{seasonDoys(1:366)}]);
    stdAsleep = nanstd([sqs_asleepDay{seasonDoys(1:366)}]);
    fprintf('%s: DAYasleep: %1.2f + %1.2f\n','All',meanAsleep/60,stdAsleep/60);
    
    meanAsleep = nanmean([sqs_asleepNight{seasonDoys(1:366)}]);
    stdAsleep = nanstd([sqs_asleepNight{seasonDoys(1:366)}]);
    fprintf('%s: NIGHTasleep: %1.2f + %1.2f\n\n','All',meanAsleep/60,stdAsleep/60);
    
    meanTrans = nanmean([sqs_trans{seasonDoys(1:366)}]);
    stdTrans = nanstd([sqs_trans{seasonDoys(1:366)}]);
    fprintf('%s: trans: %1.0f + %1.0f\n','All',meanTrans,stdTrans);
    
    meanTrans = nanmean([sqs_transDay{seasonDoys(1:366)}]);
    stdTrans = nanstd([sqs_transDay{seasonDoys(1:366)}]);
    fprintf('%s: DAYtrans: %1.0f + %1.0f\n','All',meanTrans,stdTrans);
    
    meanTrans = nanmean([sqs_transNight{seasonDoys(1:366)}]);
    stdTrans = nanstd([sqs_transNight{seasonDoys(1:366)}]);
    fprintf('%s: NIGHTtrans: %1.0f + %1.0f\n\n','All',meanTrans,stdTrans);

%% what best correlates with ODBA?
close all
f = fit(y_weather',x_odba','poly1');
figure;
plot(f,y_weather',x_odba');
[r,p] = corr(x_odba',y_weather','rows','complete');
xlabel('tempc');
ylabel('odba');
title(sprintf('TempC: r = %1.3f, p = %1.5e',r,p));

f = fit(y_dayLength',x_odba','poly1');
figure;
plot(f,y_dayLength',x_odba');
[r,p] = corr(x_odba',y_dayLength','rows','complete');
xlabel('daylength');
ylabel('odba');
title(sprintf('dayLength: r = %1.3f, p = %1.5e',r,p));

f = fit(y_asleep',x_odba','poly1');
figure;
plot(f,y_asleep',x_odba');
[r,p] = corr(x_odba',y_asleep','rows','complete');
xlabel('y_asleep');
ylabel('odba');
title(sprintf('y_asleep: r = %1.3f, p = %1.5e',r,p));


f = fit(y_asleep',y_nest','poly1');
figure;
plot(f,y_asleep',y_nest');
y_nest(y_nest == 0) = NaN;
[r,p] = corr(y_asleep',y_nest','rows','complete');
xlabel('y_asleep');
ylabel('y_nest');
title(sprintf('y_asleep: r = %1.3f, p = %1.5e',r,p));