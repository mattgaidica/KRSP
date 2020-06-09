doAnalysis = true;
weatherPath = '/Users/matt/Documents/Data/KRSP/HainesJunction_DailyTemps_Master.csv';
T_weather = readtable(weatherPath);
ssPath = '/Users/matt/Documents/Data/KRSP/SunriseSunset';
files = dir(fullfile(ssPath,'*.txt'));
T_ss = readtable(fullfile(ssPath,files(4).name));

filespath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
files = dir(fullfile(filespath,'*.meta.mat'));
if doAnalysis
    awake_sunrise = [];
    asleep_sunset = [];
    data_days = [];
    awake_length = [];
    outside_temp = [];
    for iFile = 1:numel(files)
        disp(files(iFile).name);
        load(fullfile(filespath,files(iFile).name));
        thisAwakeLength = hours(T_a.asleep-T_a.awake);
        % why are the large pos/neg numbers???
        k = thisAwakeLength > 0 & thisAwakeLength < 24;
        k_ids = find(k);
        subject_temps = [];
        for ii = 1:numel(k_ids)
            thisTemp = T_weather.Mean_Temp(find(T_weather.Date == datetime(year(T_a.sunrise(k_ids(ii))),...
                month(T_a.sunrise(k_ids(ii))),...
                day(T_a.sunrise(k_ids(ii))))));
            subject_temps = [subject_temps thisTemp];
        end
        if numel(subject_temps) > 1
            subject_temps = inpaint_nans(subject_temps);
        end
        outside_temp = [outside_temp;subject_temps'];
        awake_sunrise = [awake_sunrise;T_a.awake_sunrise(k)];
        asleep_sunset = [asleep_sunset;T_a.asleep_sunset(k)];
        data_days = [data_days;day(T_a.sunrise(k),'dayofyear')];
        awake_length = [awake_length;thisAwakeLength(k)];
    end
end

close all

ff(1200,400);
subplot(121);
x = awake_length;
y = outside_temp;
[r,p] = corr(x,y);
f = fit(x,y,'poly1');
plot(f,x,y,'k.');
xlabel('awake length');
ylabel('outside temp');
title({sprintf('r = %1.2f, p = %2.0e',r,p)});

subplot(122);
y = T_ss.day_length(data_days)/3600;
[r,p] = corr(x,y);
f = fit(x,y,'poly1');
plot(f,x,y,'k.');
xlabel('awake length');
ylabel('day length');
title({sprintf('r = %1.2f, p = %2.0e',r,p)});

[v,k] = sort(data_days);

fz = 14;
% close all
ff(900,400);
% scatter(data_days,awake_length,15,'filled','k');
f = fit(data_days,awake_length,'poly2');
plot(f,data_days,awake_length,'k.');
hold on;
plot(T_ss.day_length/3600,'r-','linewidth',2);

ylabel('hours');
xlabel('day');
xlim([1 numel(T_ss.day_length)]);
set(gca,'fontsize',fz);
