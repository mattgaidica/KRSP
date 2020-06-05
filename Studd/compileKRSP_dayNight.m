dataFiles = {'/Users/matt/Documents/Data/KRSP/J1_Sep1_2014.csv',...
    '/Users/matt/Documents/Data/KRSP/K6_Sep1_2014.csv',...
    '/Users/matt/Documents/Data/KRSP/H2_Sep1_2014.csv',...
    '/Users/matt/Documents/Data/KRSP/O1_Sep1_2014.csv'};
dayNightFiles = {'J1_Sep1_2014_meta.mat',...
    'K6_Sep1_2014_meta.mat',...
    'H2_Sep1_2014_meta.mat',...
    'O1_Sep1_2014_meta.mat'};

dayCount = 0;
odba_day = {};
odba_night = {};
temp_day = {};
temp_night = {};
days_length = [];
nights_length= [];
dates_used = {};
subjects = {};
for iFile = 1:numel(dataFiles)
    disp(dataFiles{iFile});
    A = readtable(dataFiles{iFile});
    load(dayNightFiles{iFile});
    dates = unique(A.date);
    for iDay = 1:numel(dates)
        if iDay > 1
            dayCount = dayCount + 1;
            subjects{dayCount} = dayNightFiles{iFile}(1:2);
            
            nightStart = ((iDay-2) * 86400) + (dayNight(iDay-1,2));
            dates_used{dayCount} = A.date(nightStart);
            nightEnd = dayNight(iDay,1) + ((iDay-1) * 86400);
            night_length = nightEnd - nightStart;
            nights_length(dayCount) = night_length;
            odba_night{dayCount} = A.odba(nightStart:nightEnd);
            temp_night{dayCount} = A.temp(nightStart:nightEnd);
            
            
            dayStart = dayNight(iDay,1) + (86400 * (iDay-1));
            dayEnd = (86400 * (iDay-1)) + dayNight(iDay,2);
            day_length = dayEnd - dayStart;
            days_length(dayCount) = day_length;
            odba_day{dayCount} = A.odba(dayStart:dayEnd);
            temp_day{dayCount} = A.temp(dayStart:dayEnd);
        end
    end
end
meta = struct; % squirrel ID, date, data location
meta.dates = dates_used;
meta.subjects = subjects;
meta.odba_day = odba_day;
meta.odba_night = odba_night;
meta.temp_day = temp_day;
meta.temp_night = temp_night;
meta.days_length = days_length;
meta.nights_length = nights_length;
