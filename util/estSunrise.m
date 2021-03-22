function [sunrise,sunset,noon] = estSunrise(thisYear)
% see: https://www.mathworks.com/matlabcentral/fileexchange/55509-sunrise-sunset
% noon modified from sunrise.m

lat = 60.7544541;
lon = -137.51178;
today = datetime(thisYear,1,1,'TimeZone','UTC');
useStandardTime = false;
tz = locationToTimeZone(lat,lon,useStandardTime);

yearStart = dateshift(today,"start","year");               % midnight beginning of year
yearEnd = dateshift(today,"start","year","next");          % midnight end of year
date = (yearStart+hours(12)):caldays(1):yearEnd;           % noon on each day of year

B = 360*(day(date,'dayofyear')-81)/365;
eot = minutes(9.87*sind(2*B) - 7.53*cosd(B) - 1.5*sind(B));

yearFrac = (date - yearStart) ./ (yearEnd - yearStart);    % year fraction at noon each day
gamma = 2*pi * yearFrac;                                   % year fraction in radians
delta = 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) - 0.006758*cos(2*gamma) ...
    + 0.000907*sin(2*gamma) - 0.002697*cos(3*gamma) + 0.00148*sin(3*gamma);
omega = acosd((cosd(90.833)./(cosd(lat).*cos(delta))) - tand(lat).*tan(delta));

sunrise = date - minutes(4*(lon + omega)) - eot;
sunset  = date - minutes(4*(lon - omega)) - eot;
noon = omeganoon(lon,date);

sunrise.TimeZone = tz;
sunset.TimeZone = tz;
noon.TimeZone = tz;

% % figure;
% % plot(date,timeofday(sunrise),date,timeofday(sunset),date,timeofday(noon),"LineWidth",2);
% % title("Sunrise and Sunset")
% % xlabel("Day of Year"); ylabel("Time of Day")
% % legend({'sunrise','sunset','solar noon'});
% % ylim(hours([0 24]));
end

function tz = locationToTimeZone(lat,lon,useStandardTime)
load timeZones.mat timeZones

lat = single(lat);
lon = single(lon);
for i = 1:size(timeZones,1)
    if inpolygon(lon, lat, timeZones(i).Lon,timeZones(i).Lat)
        if useStandardTime
            tz = timeZones(i).FixedID;
        else
            tz = timeZones(i).ID;
        end
        return
    end
end

tz = compose("Etc/GMT%+d",-floor((lon+7.5)/15)); % opposite sign conventions
end

function noon = omeganoon(lon,date)
% main function that computes daylength and noon time
% https://en.wikipedia.org/wiki/Sunrise_equation

% number of days since Jan 1st, 2000 12:00 UT
dte = floor(datenum(date));
n2000 = dte - datenum(2000,1,1,12,0,0) + 68.184/86400;

% mean solar moon
Js = n2000 - lon/360;

% solar mean anomaly
M = mod(357.5291 + 0.98560028*Js,360);

% center
C = 1.9148*sind(M) + 0.0200*sind(2*M) + 0.0003*sind(3*M);

% ecliptic longitude
lambda = mod(M + C + 180 + 102.9372,360);

% solar transit
Jt = 2451545.5 + Js + 0.0053*sind(M) - 0.0069*sind(2*lambda);
% % % % 
% % % % % Sun declination
% % % % delta = asind(sind(lambda)*sind(23.44));
% % % % 
% % % % % hour angle (day expressed in geometric degrees)
% % % % h = (sind(-0.83 - 2.076*sqrt(alt)/60) - sind(lat).*sind(delta))./(cosd(lat).*cosd(delta));
% % % % omega = acosd(h);
% % % % % to avoid meaningless complex angles: forces omega to 0 or 12h
% % % % omega(h<-1) = 180;
% % % % omega(h>1) = 0;
% % % % omega = real(omega);

% trick is to convert UTC here, use lat/lon afterwards
noon = datetime(Jt + datenum(2000,1,1,12,0,0) - 2451545,'ConvertFrom','datenum','TimeZone','UTC');
end
