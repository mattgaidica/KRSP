nSmooth = 500;
close all
ff(1200,600);
for iDay = 1:numel(meta.odba_day)
    subplot(211);
    plot(smooth(meta.odba_day{iDay},nSmooth));
    hold on;
    xtickslabels = xticks/3600;
    title('day odba');
    xlabel('time (h)');
    ylabel('odba');
    
    subplot(212);
    plot(smooth(meta.odba_night{iDay},nSmooth));
    hold on;
    xtickslabels = xticks/3600;
    title('night odba');
    xlabel('time (h)');
    ylabel('odba');
end