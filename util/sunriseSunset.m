function [sunrise,sunset,day_length,Tss] = sunriseSunset(t_datetime)
if ismac
    ssPath = '/Users/matt/Documents/Data/KRSP/SunriseSunset';
else
    ssPath = 'C:\Users\mgaidica\Documents\Data\KRSP\SunriseSunset';
end
Tss = readtable(fullfile(ssPath,['ss_',num2str(year(t_datetime)),'.txt']));
sunrise = Tss.sunrise(day(Tss.sunrise,'dayofyear') == day(t_datetime,'dayofyear'));
sunset = Tss.sunset(day(Tss.sunset,'dayofyear') == day(t_datetime,'dayofyear'));
day_length = Tss.day_length(day(Tss.sunrise,'dayofyear') == day(t_datetime,'dayofyear'));