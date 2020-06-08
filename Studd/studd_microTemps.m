% load('J1_Sep1_2014_meta.mat');
dates = unique(A.date);

% align data to nightStart, absolute time, stretched time?
% 12 hr = 2pi
t = linspace(0,2*pi,3600*12);
all_polar = NaN(numel(dates)-1,numel(t));
close all
h1 = ff(800,800);
h2 = ff(1200,300);
nSmooth = 100;
rlimVals = [15 29];
highAxyTemps = [];
lowAxyTemps = [];
for iDay = 1:numel(dates)
    n = find(A.date == dates(iDay));
    if iDay > 1 % night
        nightStart = ((iDay-2) * 86400) + (dayNight(iDay-1,2));
        nightEnd = dayNight(iDay,1) + ((iDay-1) * 86400);
        dayStart = dayNight(iDay,1) + (86400 * (iDay-1));
        dayEnd = (86400 * (iDay-1)) + dayNight(iDay,2);
        axyData = sum(A.odba(dayStart:dayEnd))/ day_length(iDay);
        tempData = smoothTemp(A.temp(nightStart:nightEnd),[22.01,22.5]);
%         if mean(tempData) < 24
%             continue;
%         end
        %         all_polar(iDay-1,1:numel(odbaData)) = smooth(odbaData,300);
        all_polar(iDay-1,:) = equalVectors(tempData,t);
        
        if axyData < .9
            lowAxyTemps = [lowAxyTemps;all_polar(iDay-1,:)];
            axyColor = repmat(0.5,[1,4]);
        else
            highAxyTemps = [highAxyTemps;all_polar(iDay-1,:)];
            axyColor = [1 0 0 0.5];
        end
        figure(h1);
        polarplot(t,smooth(all_polar(iDay-1,:),nSmooth),'color',axyColor);
        hold on;
        rlim(rlimVals);
        rticks(rlimVals);
        rticklabels({[num2str(rlimVals(1)),char(176),'C',],[num2str(rlimVals(2)),char(176)],'C'});
        pax = gca;
        pax.Color = [1 1 1];
        pax = gca;
        pax.ThetaZeroLocation = 'top';
        pax.ThetaDir = 'clockwise';
        pax.FontSize = 18;
        pax.Layer = 'top';
        thetaticklabels({'night start/end'});
        
        figure(h2);
        plot(smooth(all_polar(iDay-1,:),nSmooth),'color',axyColor,...
            'linewidth',1.25);
        hold on;
    end
end
% polarplot(t,normalize(nanmean(all_polar))*max(rlim),'k-','linewidth',3);

% plot(nanmean(all_polar),'k','linewidth',3);
plot(min(highAxyTemps),'r','linewidth',3);
plot(min(lowAxyTemps),'k','linewidth',3);
xlim([1 size(all_polar,2)])
ylim([17 28]);
title('temp over night (interpolated so all nights are same length)');
ylabel('temp (C)');
xlabel('time');