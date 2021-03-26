% setup with /Users/matt/Documents/MATLAB/KRSP/predict_awake.m
warning ('off','all');
weather = readtable('HainesJunction_DailyTemps_Master.csv');
warning ('on','all');
sqkey = readtable('sqkey');
filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';

Tss = makeTss(2014:2020);
unyears = unique(sqkey.year);
unyears = unyears(~isnan(unyears));
sqs_temp = NaN(numel(unyears),366);
sqs_length = NaN(numel(unyears),366);
for iYear = 1:numel(unyears)
    for iDoy = 1:366
        useId = find(weather.Year == unyears(iYear) & weather.Julian_Date == iDoy);
        if ~isempty(useId)
            sqs_temp(iYear,iDoy) = weather.Mean_Temp(useId);
        end
        useId = find(Tss.year == unyears(iYear) & Tss.doy == iDoy);
        if ~isempty(useId)
            sqs_length(iYear,iDoy) = seconds(Tss.sunset(useId)-Tss.sunrise(useId))/3600;
        end
    end
end

clc
fprintf('Mean temp: %1.2f ± %1.2f C\n',mean(nanmean(sqs_temp)),std(nanmean(sqs_temp)));

[v,k] = max(nanmean(sqs_temp));
thisDay = datestr(Tss.noon(k),'mmm dd');
fprintf('Hottest mean day: %s, doy=%i, %1.2f ± %1.2f C\n',thisDay,k,v,nanstd(sqs_temp(:,k)));

[v,k] = min(nanmean(sqs_temp));
thisDay = datestr(Tss.noon(k),'mmm dd');
fprintf('Coldest mean day: %s, doy=%i, %1.2f ± %1.2f C\n',thisDay,k,v,nanstd(sqs_temp(:,k)));

[v,k] = max(nanmean(sqs_length));
fprintf('Max Day Length: %s doy=%i, %1.2f hrs\n',datestr(Tss.noon(k),'mmm dd'),k,v);

[v,k] = min(nanmean(sqs_length));
fprintf('Min Day Length: %s doy=%i, %1.2f hrs\n',datestr(Tss.noon(k),'mmm dd'),k,v);

%% sunlight/temp polar plot (keep activity?)
% % % % % figure;
% % % % % plot(normalize(nanmean(sqs_temp),'scale'));
% % % % % hold on;
% % % % % plot(normalize(Tss.day_length,'scale'));
% % % % % plot(normalize(cellfun(@mean,sqs_odba),'scale'));

% [sunYear,sunYearAvg,sunMonth,sunHeader,monthNames] = getSunData(1:12);

sunYear = Tss.day_length(Tss.year==2016);

close all;
fh = ff(400,400);
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
edges = linspace(-pi,pi,numel(sunYear)+1);
counts = sunYear / 60 / 60 / 24;
h = polarhistogram('BinEdges',edges,'BinCounts',counts,...
    'FaceAlpha',1,'FaceColor',sunColor,'EdgeColor','none');
h = polarplot(edges,ones(size(edges)));
h.Color = 'k';

% temperature
yearlyWeather = nanmean(sqs_temp);
counts = smoothdata(normalize(yearlyWeather,'range')+1,'gaussian',100);
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
yearOdba = NaN(366,1);
for iDoy = 1:366
    useIds = find(sq_doys == iDoy);
    if numel(useIds) > 1
        yearOdba(iDoy) = mean(median(abs(sq_odba(useIds,:))));
    end
end
yearOdba = fillmissing(yearOdba,'movmean',60);
yearOdba_filt = imgaussfilt(yearOdba,10,'padding','circular');
counts = normalize(yearOdba_filt,'range')+1.5;
edges = linspace(-pi,pi,numel(counts)+1);
h = polarhistogram('BinEdges',edges,'BinCounts',counts,...
    'FaceColor','k','LineWidth',2,'FaceAlpha',0.5);
h.DisplayStyle = 'stairs';
h.EdgeColor = 'k';

% season outlines, winter begins Oct 14
seasons = linspace(-pi,pi,5)-1.8*((2*pi)/12);
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
% text(pi/2,1.1,'24 hrs','FontSize',fs,'Color','k','HorizontalAlignment','right');
text(-2.65,-1.8,'0°C','FontSize',fs,'Color','k','color','k');
text(-1.1,-2.4,'activity','FontSize',fs,'Color','k');
thetaticklabels(circshift(monthNames,6));

if true
    saveas(fh,fullfile(exportPath,'environmentalPolarPlot.png'));
    close(fh);
end
%% top table stats, see also /Users/matt/Documents/MATLAB/KRSP/Figures/cosinorEst.m
% setup seasons in /Users/matt/Documents/MATLAB/KRSP/predict_awake.m
clc

rowNames = {'Day Length (hrs)';'Sleep per day (hrs)';'Sleep in daylight (hrs)';'Sleep in darkness (hrs)';...
    'Total sleep transitions';'Sleep transitions in daylight';'Sleep transitions in darkness'};
varNames = {'Winter (Nov–Jan)','Spring (Feb–Apr)','Summer (May–Jul)','Autumn (Aug–Oct)','All'};

dayLength = {};
sleepPerDay = {};
sleepDaylight = {};
sleepDarkness = {};
sleepTrans = {};
sleepTransDaylight = {};
sleepTransDarkness = {};

for iS = 1:4
    meanDayLength = mean(Tss.day_length(seasonDoys(sIds(iS):sIds(iS+1))));
    stdDayLength = std(Tss.day_length(seasonDoys(sIds(iS):sIds(iS+1))));
    fprintf('%s: day length: %1.2f ± %1.2f\n',sTitles{iS},meanDayLength/3600,stdDayLength/3600);
    dayLength{iS} = sprintf('%1.2f ± %1.2f',meanDayLength/3600,stdDayLength/3600);
    
    meanAsleep = nanmean([sqs_asleep{seasonDoys(sIds(iS):sIds(iS+1))}]);
    stdAsleep = nanstd([sqs_asleep{seasonDoys(sIds(iS):sIds(iS+1))}]);
    fprintf('%s: asleep: %1.2f ± %1.2f\n',sTitles{iS},meanAsleep/60,stdAsleep/60);
    sleepPerDay{iS} = sprintf('%1.2f ± %1.2f',meanAsleep/60,stdAsleep/60);
    
    meanAsleep = nanmean([sqs_asleepDay{seasonDoys(sIds(iS):sIds(iS+1))}]);
    stdAsleep = nanstd([sqs_asleepDay{seasonDoys(sIds(iS):sIds(iS+1))}]);
    fprintf('%s: DAYasleep: %1.2f ± %1.2f\n',sTitles{iS},meanAsleep/60,stdAsleep/60);
    sleepDaylight{iS} = sprintf('%1.2f ± %1.2f',meanAsleep/60,stdAsleep/60);
    
    meanAsleep = nanmean([sqs_asleepNight{seasonDoys(sIds(iS):sIds(iS+1))}]);
    stdAsleep = nanstd([sqs_asleepNight{seasonDoys(sIds(iS):sIds(iS+1))}]);
    fprintf('%s: NIGHTasleep: %1.2f ± %1.2f\n\n',sTitles{iS},meanAsleep/60,stdAsleep/60);
    sleepDarkness{iS} = sprintf('%1.2f ± %1.2f',meanAsleep/60,stdAsleep/60);
    
    meanTrans = nanmean([sqs_trans{seasonDoys(sIds(iS):sIds(iS+1))}]);
    stdTrans = nanstd([sqs_trans{seasonDoys(sIds(iS):sIds(iS+1))}]);
    fprintf('%s: trans: %1.0f ± %1.0f\n',sTitles{iS},meanTrans,stdTrans);
    sleepTrans{iS} = sprintf('%1.0f ± %1.0f',meanTrans,stdTrans);
    
    meanTrans = nanmean([sqs_transDay{seasonDoys(sIds(iS):sIds(iS+1))}]);
    stdTrans = nanstd([sqs_transDay{seasonDoys(sIds(iS):sIds(iS+1))}]);
    fprintf('%s: DAYtrans: %1.0f ± %1.0f\n',sTitles{iS},meanTrans,stdTrans);
    sleepTransDaylight{iS} = sprintf('%1.0f ± %1.0f',meanTrans,stdTrans);
    
    meanTrans = nanmean([sqs_transNight{seasonDoys(sIds(iS):sIds(iS+1))}]);
    stdTrans = nanstd([sqs_transNight{seasonDoys(sIds(iS):sIds(iS+1))}]);
    fprintf('%s: NIGHTtrans: %1.0f ± %1.0f\n\n',sTitles{iS},meanTrans,stdTrans);
    sleepTransDarkness{iS} = sprintf('%1.0f ± %1.0f\n\n',meanTrans,stdTrans);
end
iS = iS + 1;
meanDayLength = mean(Tss.day_length(seasonDoys(1:366)));
stdDayLength = std(Tss.day_length(seasonDoys(1:366)));
fprintf('%s: day length: %1.2f ± %1.2f\n','All',meanDayLength/3600,stdDayLength/3600);
dayLength{iS} = sprintf('%1.2f ± %1.2f',meanDayLength/3600,stdDayLength/3600);
    
meanAsleep = nanmean([sqs_asleep{seasonDoys(1:366)}]);
stdAsleep = nanstd([sqs_asleep{seasonDoys(1:366)}]);
fprintf('%s: asleep: %1.2f ± %1.2f\n','All',meanAsleep/60,stdAsleep/60);
sleepPerDay{iS} = sprintf('%1.2f ± %1.2f',meanAsleep/60,stdAsleep/60);

meanAsleep = nanmean([sqs_asleepDay{seasonDoys(1:366)}]);
stdAsleep = nanstd([sqs_asleepDay{seasonDoys(1:366)}]);
fprintf('%s: DAYasleep: %1.2f ± %1.2f\n','All',meanAsleep/60,stdAsleep/60);
sleepDaylight{iS} = sprintf('%1.2f ± %1.2f',meanAsleep/60,stdAsleep/60);

meanAsleep = nanmean([sqs_asleepNight{seasonDoys(1:366)}]);
stdAsleep = nanstd([sqs_asleepNight{seasonDoys(1:366)}]);
fprintf('%s: NIGHTasleep: %1.2f ± %1.2f\n\n','All',meanAsleep/60,stdAsleep/60);
sleepDarkness{iS} = sprintf('%1.2f ± %1.2f',meanAsleep/60,stdAsleep/60);

meanTrans = nanmean([sqs_trans{seasonDoys(1:366)}]);
stdTrans = nanstd([sqs_trans{seasonDoys(1:366)}]);
fprintf('%s: trans: %1.0f ± %1.0f\n','All',meanTrans,stdTrans);
sleepTrans{iS} = sprintf('%1.0f ± %1.0f',meanTrans,stdTrans);

meanTrans = nanmean([sqs_transDay{seasonDoys(1:366)}]);
stdTrans = nanstd([sqs_transDay{seasonDoys(1:366)}]);
fprintf('%s: DAYtrans: %1.0f ± %1.0f\n','All',meanTrans,stdTrans);
sleepTransDaylight{iS} = sprintf('%1.0f ± %1.0f',meanTrans,stdTrans);

meanTrans = nanmean([sqs_transNight{seasonDoys(1:366)}]);
stdTrans = nanstd([sqs_transNight{seasonDoys(1:366)}]);
fprintf('%s: NIGHTtrans: %1.0f ± %1.0f\n\n','All',meanTrans,stdTrans);
sleepTransDarkness{iS} = sprintf('%1.0f ± %1.0f',meanTrans,stdTrans);

Tstats = table(dayLength',sleepPerDay',sleepDaylight',sleepDarkness',sleepTrans',sleepTransDaylight',sleepTransDarkness',...
    'VariableNames',rowNames,'RowNames',varNames);
Tstats = rows2vars(Tstats);
writetable(Tstats,'Tstats.xlsx');

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