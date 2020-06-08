
% load('KRSP_meta.mat')
% climate = readtable('/Users/matt/Documents/Data/KRSP/climate-daily.csv');
% ss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_20140101-20141231.txt');
subjects = unique(meta.subjects);
doys = [];
ys = [];
colors = lines(4);
sc = [];
day_length = [];
temps = [];
climateDays = climate.MEAN_TEMPERATURE(climate.LOCAL_YEAR == 2014);
for iDay = 1:numel(meta.dates)
    iSubject = find(strcmp(meta.subjects{iDay},subjects));
    doys(iDay) = day(meta.dates{iDay},'dayofyear');
    ys(iDay) = meta.days_length(iDay);
    sc(iDay,:) = colors(iSubject,:);
    day_length(iDay) = ss{doys(iDay),4};
    temps(iDay) = climateDays(doys(iDay));
end
% close all
ff
scatter(temps,ys,35,sc,'filled');
% scatter(temps,ys,35,sc,'filled');
xlabel('temps');
ylabel('day l');

all_years  = 2013:2019;
close all
ff;
for iYear = 1:numel(all_years)
%     plot(T.TOTAL_PRECIPITATION(T.LOCAL_YEAR == all_years(iYear)));
    plot(T.MEAN_TEMPERATURE(T.LOCAL_YEAR == all_years(iYear)));
    hold on;
end

% for iDay = 1:numel(meta.dates)
%     iSubject = find(strcmp(meta.subjects{iDay},subjects));
% end