useYears = 2014:2019;
colors = ['r','k','k','k','k','r'];
close all;
ff(1200,500);
maxs = [];
for ii = 1:numel(useYears)
    useIds = find(weather.Year == useYears(ii));
    plot(smoothdata(weather.Max_Temp(useIds),'movmean',30),'color',colors(ii),'linewidth',2);
    [v,k] = max(smoothdata(weather.Max_Temp(useIds),'movmean',30));
    maxs(ii) = k;
    hold on;
end
for ii = [2,4]
    plot([useDoys{ii}(1),useDoys{ii}(1)],ylim,'k:');
    plot([useDoys{ii}(end),useDoys{ii}(end)],ylim,'k:');
end
xlim([1 366]);
xlabel('Julian Date');
ylabel('Temp (C)');
title('Mast and Non-mast Mean Temperatures');

yyaxis right
plot(Tss.day_length,'b-');
ylabel('daylength (s)');
set(gca,'ycolor','b');

legend(compose('%i',useYears));

[v,k] = max(Tss.day_length(1:366));
fprintf("mean temp doy: %i, light doy: %i, diff: %i days\n",round(mean(maxs)),k,round(mean(maxs))-k);