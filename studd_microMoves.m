% load('J1_Sep1_2014_meta.mat');
dates = unique(A.date);

% align data to nightStart, absolute time, stretched time?
% 12 hr = 2pi
t = linspace(0,2*pi,3600*12);
all_polar = NaN(numel(dates)-1,numel(t));
close all
ff(600,600);
for iDay = 1:numel(dates)
    n = find(A.date == dates(iDay));
    if iDay > 1 % night
        nightStart = ((iDay-2) * 86400) + (dayNight(iDay-1,2));
        nightEnd = dayNight(iDay,1) + ((iDay-1) * 86400);
        odbaData = A.odba(nightStart:nightEnd);
%         all_polar(iDay-1,1:numel(odbaData)) = smooth(odbaData,300);
        all_polar(iDay-1,:) = equalVectors(odbaData,t);
        polarplot(t,all_polar(iDay-1,:));
        hold on;
%         rlim([0 rlimVal]);
        rticks([]);
        pax = gca;
        pax.Color = [1 1 1];
        pax = gca;
        pax.ThetaZeroLocation = 'top';
        pax.ThetaDir = 'clockwise';
        pax.FontSize = 18;
        pax.Layer = 'top';
        thetaticklabels(0:11);
    end
end
% polarplot(t,normalize(nanmean(all_polar))*max(rlim),'k-','linewidth',3);
ff(1200,300);
plot(smooth(nanmean(all_polar),100),'k');
xlim([1 size(all_polar,2)])
title('odba over night (sinterpolated so all nights are same length)');
ylabel('odba');
xlabel('time');