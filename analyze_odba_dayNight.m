filespath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
files = dir(fullfile(filespath,'*.meta.mat'));

doAnalysis = true;

if doAnalysis
    awake_sunrise = [];
    asleep_sunset = [];
    data_days = [];
    for iFile = 1:numel(files)
        disp(files(iFile).name);
        load(fullfile(filespath,files(iFile).name));
        awake_sunrise = [awake_sunrise;T_a.awake_sunrise];
        asleep_sunset = [asleep_sunset;T_a.asleep_sunset];
        data_days = [data_days;day(T_a.sunrise,'dayofyear')];
    end
end

binEdges = linspace(-6,6,50);
fs = 14;
close all;
ff(1200,800);
subplot(221);
histogram(awake_sunrise/3600,binEdges,'facecolor','k');
title('awake');
xlabel('before sunrise \leftarrow time (h) \rightarrow after sunrise');
ylabel('observations');
set(gca,'fontsize',fs);

subplot(223);
histogram(asleep_sunset/3600,binEdges,'facecolor','k');
title('asleep');
xlabel('before sunset \leftarrow time (h) \rightarrow after sunset');
ylabel('observations');
set(gca,'fontsize',fs);

subplot(222);
histogram(awake_sunrise(data_days < 50 | data_days > 325)/3600,binEdges);
hold on;
histogram(awake_sunrise(data_days > 100 & data_days < 260)/3600,binEdges);
title('awake');
legend({'winter','summer'});
xlabel('before sunrise \leftarrow time (h) \rightarrow after sunrise');
ylabel('observations');
set(gca,'fontsize',fs);

subplot(224);
histogram(asleep_sunset(data_days < 50 | data_days > 325)/3600,binEdges);
hold on;
histogram(asleep_sunset(data_days > 100 & data_days < 260)/3600,binEdges);
title('asleep');
legend({'winter','summer'});
xlabel('before sunset \leftarrow time (h) \rightarrow after sunset');
ylabel('observations');
set(gca,'fontsize',fs);