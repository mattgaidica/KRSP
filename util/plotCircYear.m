function plotCircYear(doys,vals,temps)

ssPath = '/Users/matt/Documents/Data/KRSP/SunriseSunset';
files = dir(fullfile(ssPath,'*.txt'));
T_ss = readtable(fullfile(ssPath,files(1).name));

useDays = [1:59 61:366]; % rm leap year?
mean_vals = NaN(365,1);
std_vals = NaN(365,1);
for ii = 1:365
    mean_vals(ii) = mean(vals(doys == useDays(ii)));
    std_vals(ii) = std(vals(doys == useDays(ii)));
end

monthNames = {'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};

op = 0.2;
fz = 14;
close all
h = ff(1200,800);
subplot(221);
polarplot(linspace(0,2*pi,365),mean_vals,'linewidth',2);
hold on;
polarplot(linspace(0,2*pi,365),mean_vals+std_vals,'linewidth',0.5,'color',[0 0 0 op]);
polarplot(linspace(0,2*pi,365),mean_vals-std_vals,'linewidth',0.5,'color',[0 0 0 op]);
polarplot(linspace(0,2*pi,numel(T_ss.day_length)),T_ss.day_length/3600,'linewidth',0.5,'color','r');
title('awake length and daylight length');

pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
pax.FontSize = fz;
pax.Layer = 'top';
rlim([0 20]);
rticks([6 max(rlim)]);
pax.Color = [1 1 1];
thetaticklabels(monthNames);

subplot(222);
histogram(vals,'facecolor','k');
ylabel('observations');
xlabel('awake length (h)');
title('awake length distribution');
set(gca,'fontsize',fz);

subplot(223);
x = T_ss.day_length(~isnan(mean_vals));
y = mean_vals(~isnan(mean_vals));
[r,p] = corr(x,y);
f = fit(x,y,'poly1');
plot(f,x,y,'k.');
title({'mean awake length vs. day length',...
    sprintf('r = %1.2f, p = %2.0e',r,p)});