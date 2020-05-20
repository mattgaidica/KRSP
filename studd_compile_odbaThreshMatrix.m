nSmooth = 45;
pThresh = 0.75;
nBins = 501;
tw = 6;

% create thresholds based on all nights for subject
if true
    subjects = unique(meta.subjects);
    compiled_night_odba = cell(numel(subjects),1);
    compiled_day_odba = cell(numel(subjects),1);
    for iDay = 1:numel(meta.dates)
        iSubject = find(strcmp(meta.subjects{iDay},subjects));
        odba_night = meta.odba_night{iDay};
        compiled_night_odba{iSubject} = [compiled_night_odba{iSubject};odba_night];
        odba_day = meta.odba_day{iDay};
        compiled_day_odba{iSubject} = [compiled_day_odba{iSubject};odba_day];
    end
    v_nThresh = [];
    v_dThresh = [];
    v_dMean = [];
    v_nMean = [];
    v_dStd = [];
    v_nStd = [];
    for iSubject = 1:numel(subjects)
        v_night = sort(compiled_night_odba{iSubject});
        v_nThresh(iSubject) = v_night(round(numel(v_night)*pThresh));
        v_day = sort(compiled_day_odba{iSubject});
        v_dThresh(iSubject) = v_day(round(numel(v_day)*pThresh));
        v_dMean(iSubject) = mean(v_day);
        v_nMean(iSubject) = mean(v_night);
        v_dStd(iSubject) = std(v_day);
        v_nStd(iSubject) = std(v_night);
    end
end

close all
ff(1200,900);
rows = 6;
cols = 2;
caxisVals = [0 1.5];
nTicks = 15;

% binEdges = linspace(0,-tw*3600,nBins);
odba_hist = [];
for iDay = 1:numel(meta.dates)
    iSubject = find(strcmp(meta.subjects{iDay},subjects));
    odba_day = smooth(meta.odba_day{iDay}(1:tw*3600),nSmooth);
%     odba_ids = find(odba_night(end+tWake*3600:end) > v_nThresh(iSubject));
    odba_hist(iDay,:) = equalVectors((odba_day-v_dMean(iSubject))/v_dStd(iSubject),binEdges(1:end-1));%histcounts(odba_ids,binEdges);
end
[~,k] = sort(sum(odba_hist'));
subplot(rows,cols,[1,3]);
imagesc(odba_hist(k,:));
colormap(hot);
xticks(linspace(1,size(odba_hist,2),nTicks));
xticklabels(compose('%1.2f',linspace(0,tw,nTicks)));
xtickangle(60);
xlabel('time rel. to wake (h)');
yticks([1,size(odba_hist,1)]);
ylabel('day (sorted by odba)');
caxis(caxisVals);
xlabel('time from start of day');
title('day');

subplot(rows,cols,5);
plot(sum(odba_hist),'k','linewidth',2);
xlim([1,size(odba_hist,2)]);
xticks(linspace(1,size(odba_hist,2),nTicks));
xticklabels(compose('%1.2f',linspace(0,tw,nTicks)));
xtickangle(60);
ylabel('odba z');
xlabel('time from start of day');



odba_hist = [];
for iDay = 1:numel(meta.dates)
    iSubject = find(strcmp(meta.subjects{iDay},subjects));
    odba_day = smooth(meta.odba_day{iDay}(end-tw*3600:end),nSmooth);
%     odba_ids = find(odba_night(end+tWake*3600:end) > v_nThresh(iSubject));
    odba_hist(iDay,:) = equalVectors((odba_day-v_dMean(iSubject))/v_dStd(iSubject),binEdges(1:end-1));%histcounts(odba_ids,binEdges);
end
[~,k] = sort(sum(odba_hist'));
subplot(rows,cols,[2,4]);
imagesc(odba_hist(k,:));
colormap(hot);
nTicks = 15;
xticks(linspace(1,size(odba_hist,2),nTicks));
xticklabels(compose('%1.2f',linspace(-tw,0,nTicks)));
xtickangle(60);
xlabel('time rel. to wake (h)');
yticks([1,size(odba_hist,1)]);
ylabel('day (sorted by odba)');
caxis(caxisVals);
xlabel('time from end of day');
title('day');

subplot(rows,cols,6);
plot(sum(odba_hist),'k','linewidth',2);
xlim([1,size(odba_hist,2)]);
xticks(linspace(1,size(odba_hist,2),nTicks));
xticklabels(compose('%1.2f',linspace(-tw,0,nTicks)));
xtickangle(60);
ylabel('odba z');
xlabel('time from end of day');



odba_hist = [];
for iDay = 1:numel(meta.dates)
    iSubject = find(strcmp(meta.subjects{iDay},subjects));
    odba_night = smooth(meta.odba_night{iDay}(1:tw*3600),nSmooth);
%     odba_ids = find(odba_night(end+tWake*3600:end) > v_nThresh(iSubject));
    odba_hist(iDay,:) = equalVectors((odba_night-v_nMean(iSubject))/v_nStd(iSubject),binEdges(1:end-1));%histcounts(odba_ids,binEdges);
end
[~,k] = sort(sum(odba_hist'));
subplot(rows,cols,[7,9]);
imagesc(odba_hist(k,:));
colormap(hot);
nTicks = 15;
xticks(linspace(1,size(odba_hist,2),nTicks));
xticklabels(compose('%1.2f',linspace(0,tw,nTicks)));
xtickangle(60);
xlabel('time rel. to wake (h)');
yticks([1,size(odba_hist,1)]);
ylabel('day (sorted by odba)');
caxis(caxisVals);
xlabel('time from start of night');
title('night');

subplot(rows,cols,11);
plot(sum(odba_hist),'k','linewidth',2);
xlim([1,size(odba_hist,2)]);
xticks(linspace(1,size(odba_hist,2),nTicks));
xticklabels(compose('%1.2f',linspace(0,tw,nTicks)));
xtickangle(60);
ylabel('odba count');
xlabel('time from start of night');



odba_hist = [];
for iDay = 1:numel(meta.dates)
    iSubject = find(strcmp(meta.subjects{iDay},subjects));
    odba_night = smooth(meta.odba_night{iDay}(end-tw*3600:end),nSmooth);
%     odba_ids = find(odba_night(end+tWake*3600:end) > v_nThresh(iSubject));
    odba_hist(iDay,:) = equalVectors((odba_night-v_nMean(iSubject))/v_nStd(iSubject),binEdges(1:end-1));%histcounts(odba_ids,binEdges);
end
[~,k] = sort(sum(odba_hist'));
subplot(rows,cols,[8,10]);
imagesc(odba_hist(k,:));
colormap(hot);
nTicks = 15;
xticks(linspace(1,size(odba_hist,2),nTicks));
xticklabels(compose('%1.2f',linspace(-tw,0,nTicks)));
xtickangle(60);
xlabel('time rel. to wake (h)');
yticks([1,size(odba_hist,1)]);
ylabel('day (sorted by odba)');
caxis(caxisVals);
xlabel('time from end of night');
title('night');

subplot(rows,cols,12);
plot(sum(odba_hist),'k','linewidth',2);
xlim([1,size(odba_hist,2)]);
xticks(linspace(1,size(odba_hist,2),nTicks));
xticklabels(compose('%1.2f',linspace(-tw,0,nTicks)));
xtickangle(60);
ylabel('odba count');
xlabel('time from end of night');