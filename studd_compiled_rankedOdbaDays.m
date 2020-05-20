pThresh = 0.75;
nSmooth = 45;
minpeakdist = 45;
% create thresholds based on all nights for subject
if pThresh ~= saveThresh
    disp('redoing thesholds...');
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
    for iSubject = 1:numel(subjects)
        v_night = sort(compiled_night_odba{iSubject});
        v_nThresh(iSubject) = v_night(round(numel(v_night)*pThresh));
        v_day = sort(compiled_day_odba{iSubject});
        v_dMean(iSubject) = mean(v_day);
        v_dThresh(iSubject) = v_day(round(numel(v_day)*pThresh));
    end
    saveThresh = pThresh; % init saveThresh manually
end

% close all
ff(1200,600);
odba_day_sum = [];
odba_night_sum = [];
for iDay = 1:numel(meta.dates)
    iSubject = find(strcmp(meta.subjects{iDay},subjects));
    
    odba_day = smooth(meta.odba_day{iDay},nSmooth);
    odba_day_sum(iDay) = sum(odba_day - v_dMean(iSubject)) / meta.days_length(iDay);
%     [locs, pks] = peakseek(odba_day,minpeakdist,v_dThresh(iSubject));
%     odba_day_sum(iDay) = numel(locs) / meta.days_length(iDay);
    if iDay == 1
        subplot(231);
        plot(odba_day);
        hold on;
        plot([min(xlim),max(xlim)],[v_dMean(iSubject),v_dMean(iSubject)],'k--');
        title({'1 day, smoothed day ODBA w/ mean_{subject}',...
            sprintf('thresh = %1.2f, minpeakdist = %is, smooth = %is',pThresh,minpeakdist,nSmooth)});
        xticklabels(compose('%1.2f',xticks/3600));
        xtickangle(60);
        xlim([1,numel(odba_day)]);
        xlabel('time (h)');
        ylabel('odba');
    end
    
    odba_night = smooth(meta.odba_night{iDay},nSmooth); % (end-3600*3:end)
%     odba_night_sum(iDay) = sum(odba_night > v_nThresh(iSubject)) / meta.nights_length(iDay);
    [locs, pks] = peakseek(odba_night,minpeakdist,v_nThresh(iSubject));
    odba_night_sum(iDay) = numel(locs) / meta.nights_length(iDay);
    if iDay == 1
        subplot(234);
        plot(odba_night);
        hold on;
        plot(locs,pks,'r*');
        plot([min(xlim),max(xlim)],[v_nThresh(iSubject),v_nThresh(iSubject)],'k--');
        title('1 night, smoothed night ODBA awakenings');
        xticklabels(compose('%1.2f',xticks/3600));
        xtickangle(60);
        xlim([1,numel(odba_night)]);
        xlabel('time (h)');
        ylabel('odba');
    end
end

lw = 2;
subplot(232);
[vd,kd] = sort(odba_day_sum);
plot(odba_day_sum(kd),'k','linewidth',lw);
ylabel('day odba / time');
yyaxis right;
plot(odba_night_sum(kd),'linewidth',lw);
ylabel('night odba / time');
xlim([size(odba_night_sum)]);
xlabel('day');
title('odba sorted by day');

subplot(235);
[vn,kn] = sort(odba_night_sum);
plot(odba_day_sum(kn),'k','linewidth',lw);
ylabel('day odba / time');
yyaxis right;
plot(odba_night_sum(kn),'linewidth',lw);
ylabel('night odba / time');
xlim([size(odba_night_sum)]);
xlabel('day');
title('odba sorted by night');

[A,I] = rmoutliers(odba_day_sum);
B = odba_night_sum(~I);
[B,I] = rmoutliers(B);
A = A(~I);
[r,p] = corr(A',B');
subplot(233);
% f = fit(odba_day_sum',odba_night_sum','poly1');
% plot(f,A,B); hold on;
f = fit(A',B','poly1');
plot(f,A,B);
xlabel('day');
ylabel('night');
title(sprintf('r = %1.2f, p = %1.5f (n=%i->%i)',r,p,numel(odba_day_sum),numel(B)));

dayrank = [];
nightrank = [];
for iDay = 1:numel(odba_day_sum)
    dayrank(iDay) = find(odba_day_sum(iDay) == vd,1,'first');
    nightrank(iDay) = find(odba_night_sum(iDay) == vn,1,'first');
end
[r,p] = corr(dayrank',nightrank');
subplot(236);
f = fit(dayrank',nightrank','poly1');
plot(f,dayrank,nightrank);
xlabel('day rank');
ylabel('night rank');
title(sprintf('r = %1.2f, p = %1.5f',r,p));