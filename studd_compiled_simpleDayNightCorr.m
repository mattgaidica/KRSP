if true
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
    v_nMean = [];
    for iSubject = 1:numel(subjects)
        v_day = sort(compiled_day_odba{iSubject});
        v_night = sort(compiled_night_odba{iSubject});
        v_dThresh(iSubject) = v_day(round(numel(v_day)*pThreshs(iiN)));
        v_nThresh(iSubject) = v_night(round(numel(v_night)*pThreshs(iiN)));
        v_dMean(iSubject) = mean(v_day);
        v_nMean(iSubject) = mean(v_night);
    end
end
for iDay = 1:numel(meta.dates)
    iSubject = find(strcmp(meta.subjects{iDay},subjects));
    odba_day = smooth(meta.odba_day{iDay},nSmooth);
    odba_day_counts(iDay) = sum(odba_day - v_dMean(iSubject)) / meta.days_length(iDay);
    odba_night = smooth(meta.odba_night{iDay},nSmooth);
    odba_night_counts(iDay) = sum(odba_night - v_nMean(iSubject)) / meta.nights_length(iDay);
end

ff(800,300);
% w/ outliers
[r,p] = corr(odba_day_counts',odba_night_counts');
f = fit(odba_day_counts',odba_night_counts','poly1');
subplot(121);
plot(f,odba_day_counts',odba_night_counts','k.');
title(sprintf('r = %1.2f, p = %1.5f, days = %i',r,p,numel(odba_day_counts)));
xlabel('day ODBA');
ylabel('night ODBA');

% rmoutliers
[A,I] = rmoutliers(odba_day_counts);
B = odba_night_counts(~I);
[B,I] = rmoutliers(B);
A = A(~I);
[r,p] = corr(A',B');
f = fit(A',B','poly1');
subplot(122);
plot(f,A',B','k.');
title(sprintf('r = %1.2f, p = %1.5f, days = %i',r,p,numel(A)));
xlabel('day ODBA');
ylabel('night ODBA');