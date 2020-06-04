function [sunrise,sunset] = sunriseSunset(t_datetime)
if ismac
    ssPath = '/Users/matt/Documents/Data/KRSP/SunriseSunset';
else
    ssPath = 'C:\Users\mgaidica\Documents\Data\KRSP\SunriseSunset';
end
T = readtable(fullfile(ssPath,['ss_',num2str(year(t_datetime)),'.txt']));
[v,k] = min(abs(T.sunrise - t_datetime));
sunrise = T.sunrise(k);
sunset = T.sunset(k);
% day_length = T.day_length;