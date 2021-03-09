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
        if ~isempty(sqkey.filename{iSq})% && ~any(ismember(sqkey.year(iSq),[2014,2019])) % ~(strcmp(sqkey.source{iSq},'ES') &&
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            if sqkey.isValid(iSq)
                T = detect_sleepWake2(T,60);
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
    end
    do = false;
end

%%
% % % % % figure;
% % % % % plot(normalize(nanmean(sqs_temp),'scale'));
% % % % % hold on;
% % % % % plot(normalize(Tss.day_length,'scale'));
% % % % % plot(normalize(cellfun(@mean,sqs_odba),'scale'));

[sunYear,sunYearAvg,sunMonth,sunHeader,monthNames] = getSunData(1:12);

close all;
ff(400,400);
rlimVal = 2.65;
colors = lines(8);
sunColor = colors(3,:);
accelColor = colors(5,:); %colors(5,:);

% night
edges = linspace(-pi,pi,13);
counts = ones(1,12);
h = polarhistogram('BinEdges',edges,'BinCounts',counts,...
    'FaceColor','k','LineWidth',0.25,'FaceAlpha',0.9);
h.EdgeColor = 'none';
hold on;

% day
edges = linspace(-pi,pi,size(sunYear,1)+1);
counts = sunYear(:,4) / 60 / 24;
h = polarhistogram('BinEdges',edges,'BinCounts',counts,...
    'FaceAlpha',1,'FaceColor',sunColor,'EdgeColor','none');
h = polarplot(edges,ones(size(edges)));
h.Color = 'k';

% temperature
yearlyWeather = nanmean(sqs_temp);
counts = smoothdata(normalize(yearlyWeather,'range')+1,'gaussian',50);
edges = linspace(-pi,pi,numel(counts));
colors = parula(numel(counts));
colorLookup = linspace(1,2,size(colors,1));
for ii = 1:numel(counts)
    %     ct = zeros(size(counts));
    %     ct(ii) = counts(ii);
    colorId = closest(colorLookup,counts(ii));
    %     h = polarhistogram('BinEdges',edges,'BinCounts',ct,...
    %         'FaceColor',colors(colorId,:),'EdgeColor','none','FaceAlpha',0.5);
    h = polarplot([edges(ii) edges(ii)],[1 counts(ii)],'-','Color',colors(colorId,:),...
        'MarkerSize',15,'LineWidth',1);
    hold on;
end
h = polarplot(edges,ones(size(edges))*1.5,'linewidth',1);
h.Color = repmat(0,[1,3]);
h.LineStyle = ':';

% odba
yearOdba = fillmissing(cellfun(@mean,sqs_odba),'movmean',40);
yearOdba_filt = imgaussfilt(yearOdba,3,'padding','circular');
counts = normalize(yearOdba_filt,'range')+1.5;
edges = linspace(-pi,pi,numel(counts)+1);
h = polarhistogram('BinEdges',edges,'BinCounts',counts,...
    'FaceColor','k','LineWidth',2,'FaceAlpha',0.5);
h.DisplayStyle = 'stairs';
h.EdgeColor = 'k';

% season outlines
seasons = linspace(-pi,pi,5)-2*((2*pi)/12);
n = 1000;
spacing = pi/128;
for iS = 1:4
    theseTheta = linspace(seasons(iS)+spacing,seasons(iS+1)-spacing,n);
    theseRho = ones(1,n)*rlimVal;
    h = polarplot(theseTheta,theseRho,'color',repmat(0.7+(iS/17),[1,3]),'LineWidth',8);
end

fs = 14;
pax = gca;
pax.ThetaZeroLocation = 'bottom';
pax.ThetaDir = 'clockwise';
pax.FontSize = fs;
pax.Layer = 'top';
rlim([0 rlimVal]);
rticks([]);
pax.Color = [1 1 1];
% rticklabels({'','',''});
text(pi/2,1.1,'24 hrs','FontSize',fs,'Color','k','HorizontalAlignment','right');
text(-2.65,1.67,'0°C','FontSize',fs,'Color','k','color','k');
text(0,2.31,'activity','FontSize',fs,'Color','k');
thetaticklabels(circshift(monthNames,6));

%% top table stats, see also /Users/matt/Documents/MATLAB/KRSP/Figures/cosinorEst.m
clc
sTitles = {'winter','spring','summer','fall'};
sIds = round(linspace(1,366,5));
seasonDoys = circshift(1:366,60);

for iS = 1:4
    meanDayLength = mean(Tss.day_length(seasonDoys(sIds(iS):sIds(iS+1))));
    stdDayLength = std(Tss.day_length(seasonDoys(sIds(iS):sIds(iS+1))));
    fprintf('%s: day length: %1.2f + %1.2f\n',sTitles{iS},meanDayLength/3600,stdDayLength/3600);
    
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
meanDayLength = mean(Tss.day_length(seasonDoys(1:366)));
stdDayLength = std(Tss.day_length(seasonDoys(1:366)));
fprintf('%s: day length: %1.2f + %1.2f\n','All',meanDayLength/3600,stdDayLength/3600);
    
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