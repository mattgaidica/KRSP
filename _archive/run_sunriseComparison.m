ssPath = '/Users/matt/Documents/Data/KRSP/SunriseSunset';
files = dir(fullfile(ssPath,'*.txt'));
colors = lines(numel(files));
close all
ff(800,600);
fz = 14;
lns = [];
for iFile = 1:numel(files)
    disp(files(iFile).name);
    T = readtable(fullfile(ssPath,files(iFile).name));
    
    subplot(211);
    lns(iFile) = plot(T.day_length,'color',colors(iFile,:),'linewidth',2);
    hold on;
    title('day length');
    yticklabels(compose('%1.1f',yticks/3600));
    ylabel('hours');
    xlabel('day');
    set(gca,'fontsize',fz);
    xlim([1 numel(T.day_length)]);
    
    subplot(212);
    plot(timeofday(T.sunrise),'color',colors(iFile,:),'linewidth',2);
    hold on;
    plot(timeofday(T.sunset),'color',colors(iFile,:),'linewidth',2);
    title('sunrise/sunset');
    ylabel('time of day');
    xlabel('day');
    set(gca,'fontsize',fz);
    xlim([1 numel(T.day_length)]);
end
subplot(211);
legend(lns,{files.name},'location','east');
legend box off;