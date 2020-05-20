function dayNight = markDayNight(A)

dates = unique(A.date);
dayNight = zeros(numel(dates),2);
for iDay = 1:numel(dates)
    h = ff(1200,300);
    n = (A.date == dates(iDay));
    plot(smoothTemp(A.temp(n),[22.01,22.5]),'k');
    ylim([10 27]);
    hold on;
    yyaxis right;
    plot(smooth(A.odba(n),200),'r');
    title(datestr(dates(iDay)));
    ylim([0 4]);
    legend({'temp','axy'});
    [xs,~] = ginput(2);
    dayNight(iDay,:) = round(xs);
    close(h);
end