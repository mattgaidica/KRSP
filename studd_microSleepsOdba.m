nSweep = 50;
tw = 3600/2;
fromWakeup = true;
useHours = 6;
if fromWakeup
    sweep = fliplr(round(linspace(1,useHours*3600,nSweep)));
else
    sweep = round(linspace(1,useHours*3600,nSweep));
end
pThreshs = linspace(0.5,0.99,10);
% pThreshs = .75;
colors = cool(numel(pThreshs));
nSmooth = 45;
minpeakdist = 45;

% create thresholds based on all nights for subject
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
    for iSubject = 1:numel(subjects)
        v_day = sort(compiled_day_odba{iSubject});
        v_night = sort(compiled_night_odba{iSubject});
        for iiN = 1:numel(pThreshs)
            v_dThresh(iSubject,iiN) = v_day(round(numel(v_day)*pThreshs(iiN)));
            v_nThresh(iSubject,iiN) = v_night(round(numel(v_night)*pThreshs(iiN)));
        end
        v_dMean(iSubject) = mean(v_day);
    end
end

% close all
h1 = ff(1200,600);
lns = [];
for jj = 1:numel(pThreshs)
    rs = [];
    pvals = [];
    for ii = 1:nSweep
        odba_night_counts = [];
        odba_day_counts = [];
        fprintf('ii %i, jj %i\n',ii,jj);
        for iDay = 1:numel(meta.dates)
            iSubject = find(strcmp(meta.subjects{iDay},subjects));
            
            odba_day = smooth(meta.odba_day{iDay},nSmooth);
            odba_day_counts(iDay) = sum(odba_day - v_dMean(iSubject)) / meta.days_length(iDay);
            if fromWakeup
                odba_night = smooth(meta.odba_night{iDay}(end-sweep(ii)-tw:end-sweep(ii)+1),nSmooth);
                sw = -1;
            else
                odba_night = smooth(meta.odba_night{iDay}(sweep(ii):sweep(ii)+tw-1),nSmooth);
                sw = 1;
            end
            [locs, pks] = peakseek(odba_night,minpeakdist,v_nThresh(iSubject,jj));
            odba_night_counts(iDay) = numel(locs);
            %             odba_night_counts(iDay) = sum(odba_night(end-sweep(ii)-tw:end-sweep(ii)) > v_nThresh(jj));
            %             odba_counts(iDay) = sum(odba_night(sweep(ii):sweep(ii)+tw) > pThresh);
        end
        [A,I] = rmoutliers(odba_day_counts);
        B = odba_night_counts(~I);
        [B,I] = rmoutliers(B);
        A = A(~I);
        
        [r,pval] = corr(A',B');
        rs(ii) = r;
        pvals(ii) = pval;
    end
    lw = 1.5;
    opac = 0.25;
    figure(h1);
    subplot(211);
    plot(sw*sweep/3600,rs,'color',[colors(jj,:),opac],'linewidth',lw);
    hold on;
    plot(sw*sweep(find(pvals < 0.05))/3600,rs(pvals < 0.05),'*','color',[colors(jj,:),opac]);
    title(sprintf('axy corr (night-day) %1.2f hour window',tw/3600));
    ylabel('r');
    if fromWakeup
        xlabel('time from wakeup');
    else
        xlabel('time after going to sleep');
    end
    
    subplot(212);
    lns(jj) = plot(sw*sweep/3600,pvals,'color',[colors(jj,:),opac],'linewidth',lw);
    hold on;
    plot(sw*sweep(find(pvals < 0.05))/3600,pvals(pvals < 0.05),'*','color',[colors(jj,:),opac]);
    ylabel('p-value');
    if fromWakeup
        xlabel('time from wakeup');
    else
        xlabel('time after going to sleep');
    end
    drawnow;
end
legend(lns,compose('%1.2f',pThreshs),'location','west');
% legend(lns(1:8:50),compose('%1.2f',pThreshs(1:8:50)),'location','west');
legend box off;