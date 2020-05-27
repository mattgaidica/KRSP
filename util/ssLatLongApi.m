function T = ssLatLongApi(firstDate,lastDate)
% firstDate = datetime(2019,1,1);
% lastDate = datetime(2019,12,31);
% lat/lng for Haines Junction, Yukon, CA
savePath = '/Users/matt/Documents/Data/KRSP/SunriseSunset';
latStr = '60.754149';
lngStr = '-137.510827';
tz = 'America/Whitehorse';
options = weboptions('Timeout',10);
header = {'sunrise',...
    'sunset',...
    'solar_noon',...
    'day_length',...
    'civil_twilight_begin',...
    'civil_twilight_end',...
    'nautical_twilight_begin',...
    'nautical_twilight_end',...
    'astronomical_twilight_begin',...
    'astronomical_twilight_end'};
varTypes = {'datetime','datetime','datetime','double',...
    'datetime','datetime','datetime','datetime','datetime','datetime'};
nDays = daysact(firstDate,lastDate) + 1;
T = table('Size',[nDays,numel(header)],'VariableTypes',varTypes,'VariableNames',header);
filename = ['ss_',datestr(datetime(firstDate),'yyyymmdd'),'-',datestr(datetime(lastDate),'yyyymmdd')];

for iDay = 1:nDays
    url = ssUrl(latStr,lngStr,datestr(datetime(firstDate+iDay-1),'yyyy-mm-dd'));
    disp(url);
    try
        S = webread(url,options);
    catch ME
        error('Unable to complete web request.');
        S = [];
    end
    pause(0.5); % throttle

    for iField = 1:numel(header)
        result = getfield(S.results,header{iField});
        if strcmp(header{iField},'day_length')
            t = result;
        else
            t = datetime(result,'TimeZone',tz,'InputFormat','yyyy-MM-dd''T''HH:mm:ssX');
            t = datetime(t,'TimeZone',''); % unzone the entry for table
        end
        T{iDay,iField} = t;
    end
    if mod(iDay,10) == 0
        writetable(T,fullfile(savePath,filename));
        disp(['saving: ',filename]);
    end
end
writetable(T,fullfile(savePath,filename));
end

function url = ssUrl(latStr,lngStr,dayDate)
endpoint = 'https://api.sunrise-sunset.org';
url = [endpoint,'/json?lat=',latStr,'&lng=',lngStr,'&date=',dayDate,'&formatted=0'];
end