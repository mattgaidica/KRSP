% load('J1_Sep1_2014_meta.mat');
% data starts and ends with partial night
% therefore, days = nights + 1
dates = unique(A.date);

day_temp = NaN(size(dates));
day_axy = day_temp;
day_length = day_temp;
day_temp_var = day_temp;
day_axy_var = day_temp;
night_temp = day_temp;
night_axy = day_temp;
night_length = day_temp;
night_temp_var = day_temp;
night_axy_var = day_temp;
% night = "night before"
% day = "day of"
for iDay = 1:numel(dates)
    n = find(A.date == dates(iDay));
    if iDay > 1 % night
        nightStart = ((iDay-2) * 86400) + (dayNight(iDay-1,2));
        nightEnd = dayNight(iDay,1) + ((iDay-1) * 86400);
        night_length(iDay) = nightEnd - nightStart;
        night_temp(iDay) = sum(smoothTemp(A.temp(nightStart:nightEnd),[22.01,22.5])) / night_length(iDay);
        night_axy(iDay) = sum(A.odba(nightStart:nightEnd)) / night_length(iDay);
        night_temp_var(iDay) = std(A.temp(nightStart:nightEnd));
        night_axy_var(iDay) = std(A.odba(nightStart:nightEnd));
    end
    dayStart = dayNight(iDay,1) + (86400 * (iDay-1));
    dayEnd = (86400 * (iDay-1)) + dayNight(iDay,2);
    day_length(iDay) = dayEnd - dayStart;
    day_temp(iDay) = sum(smoothTemp(A.temp(dayStart:dayEnd),[22.01,22.5])) / day_length(iDay);
    day_axy(iDay) = sum(A.odba(dayStart:dayEnd))/ day_length(iDay);
    day_temp_var(iDay) = std(A.temp(dayStart:dayEnd));
    day_axy_var(iDay) = std(A.odba(dayStart:dayEnd));
end

lw = 2;
close all

ff(1200,400);
subplot(121);
[r,pval] = corr(day_axy(2:end),night_temp(2:end));
f = fit(day_axy(2:end),night_temp(2:end),'poly1');
scatter(day_axy(2:end),night_temp(2:end),35,'filled');
hold on;
plot(f);
title({'prev-night temp vs. day-of axy per-minute', sprintf('r=%1.3f, p=%1.3f',r,pval)});
xlabel('day-of axy');
ylabel('prev-night temp');
subplot(122);
[r,pval] = corr(day_axy(2:end-1),night_temp(3:end));
f = fit(day_axy(2:end-1),night_temp(3:end),'poly1');
scatter(day_axy(2:end-1),night_temp(3:end),35,'filled');
title({'night-of temp vs. day-of axy per-minute', sprintf('r=%1.3f, p=%1.3f',r,pval)});
hold on;
plot(f);
xlabel('day-of axy');
ylabel('night-of temp');


% rmoutliers code
t_night_axy = night_axy(2:end);
t_day_axy = day_axy(2:end);
[night_axy_rmoutliers,I] = rmoutliers(t_night_axy);
day_axy_rmoutliers = t_day_axy(~I);

ff(1200,400);
subplot(121);
[r,pval] = corr(day_axy_rmoutliers,night_axy_rmoutliers);
f = fit(day_axy_rmoutliers,night_axy_rmoutliers,'poly1');
scatter(day_axy_rmoutliers,night_axy_rmoutliers,35,'filled');
hold on;
plot(f);
title({'prev-night vs. day-of axy per-minute', sprintf('r=%1.3f, p=%1.3f',r,pval)});
xlabel('day-of axy');
ylabel('prev-night axy');
subplot(122);
[r,pval] = corr(day_axy(2:end-1),night_axy(3:end));
f = fit(day_axy(2:end-1),night_axy(3:end),'poly1');
scatter(day_axy(2:end-1),night_axy(3:end),35,'filled');
title({'night-of vs. day-of axy per-minute', sprintf('r=%1.3f, p=%1.3f',r,pval)});
hold on;
plot(f);
xlabel('day-of axy');
ylabel('night-of axy');

ff(1200,800);
subplot(311)
plot(day_temp,'linewidth',lw);
hold on;
plot(night_temp,'linewidth',lw);
xlim([1 numel(dates)]);
xticks(1:numel(dates));
legend({'day','night'});
legend box off;
grid on;
title('sum temp per-second');

subplot(312)
plot(day_axy,'linewidth',lw);
yyaxis right;
plot(night_axy,'linewidth',lw);
title('sum axy per-second');
xlim([1 numel(dates)]);
xticks(1:numel(dates));
grid on;
legend({'day','night'});
legend box off;

subplot(313)
plot(day_temp_var,'linewidth',lw);
hold on;
plot(night_temp_var,'linewidth',lw);
yyaxis right;
plot(day_axy_var,':','linewidth',lw);
plot(night_axy_var,'g:','linewidth',lw);
xlim([1 numel(dates)]);
xticks(1:numel(dates));
grid on;
title('temp var');
legend({'day temp','night temp','day axy','night axy'});
legend box off;


% figure;
% scatter(day_length/3600,night_length/3600,'filled',25);


ff(800,300);
plot(day_length/3600,'linewidth',lw);
hold on;
plot(night_length/3600,'linewidth',lw);

legend({'day','night'});
legend box off;
xlim([1 numel(dates)]);
xticks(1:numel(dates));
ylabel('hours');
title('day/night length');
grid on;