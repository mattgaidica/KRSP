% load('KRSP_meta.mat')
nSmooth = 45;
subjects = unique(meta.subjects);
% colors = parula(numel(meta.dates));
colors = lines(numel(subjects));
close all
ff;
op = 0.3;
for iDay = 1:numel(meta.dates)
    iSubject = find(strcmp(meta.subjects{iDay},subjects));
    
    odba_day = smooth(meta.odba_day{iDay},nSmooth);
    odba_night = smooth(meta.odba_night{iDay},nSmooth);
    plot(0:numel(odba_day)-1,cumsum(odba_day),'color',[colors(iSubject,:),op]);
    hold on
    plot(-numel(odba_night)+1:0,cumsum(flip(odba_night)) - max(cumsum(odba_night)),...
        'color',[colors(iSubject,:),op]);
end
xticks([min(xlim) 0 max(xlim)]);
yticks([min(ylim) 0 max(ylim)]);
grid on;