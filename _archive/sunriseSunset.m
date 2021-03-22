function [sunrise,sunset,day_length,Tss_all] = sunriseSunset(t_datetime,Tss_all)
if ismac
    ssPath = '/Users/matt/Documents/Data/KRSP/SunriseSunset';
else
    ssPath = 'C:\Users\mgaidica\Documents\Data\KRSP\SunriseSunset';
end

startYear = 2013;
files = dir(fullfile(ssPath,'*.txt'));
if isempty(Tss_all)
    Tss_all = {};
    for iFile = 1:numel(files)
        thisYear = str2double(files(iFile).name(4:4+3));
        Tss_all{thisYear-startYear+1} = readtable(fullfile(ssPath,files(iFile).name));
    end
end

theseYears = year(t_datetime);
theseDoys = day(t_datetime,'dayofyear');
sunrise = NaT(0);
sunset = NaT(0);
day_length = [];
for iYear = unique(theseYears)'
    sunrise = [sunrise;Tss_all{iYear-startYear+1}.sunrise(theseDoys(theseYears==iYear))];
    sunset = [sunset;Tss_all{iYear-startYear+1}.sunset(theseDoys(theseYears==iYear))];
    day_length = [day_length;Tss_all{iYear-startYear+1}.day_length(theseDoys(theseYears==iYear))];
end