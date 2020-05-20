nSmooth = 45;
showHours = 5;
nBins = 100;
allDiffs = [];
pThreshs = .75;%linspace(.73,0.77,3);
colors = [0 0 0]; %cool(numel(pThreshs));
% % % % saveSubject = meta.subjects{1};
close all
ff(600,600);
lns = [];
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
for iPt = 1:numel(pThreshs)
    for iDay = 1:numel(meta.dates)
        disp(iDay);
        iSubject = find(strcmp(meta.subjects{iDay},subjects));
        
        odba_night = smooth(meta.odba_night{iDay},nSmooth);
        [locs, pks] = peakseek(odba_night,nSmooth,v_nThresh(iSubject,iPt));
        if false
            ff(600,400);
            plot(odba_night);
            hold on;
            plot(locs,pks,'r*');
            plot([min(xlim),max(xlim)],[v_nThresh(iSubject,iPt),v_nThresh(iSubject,iPt)],'k--');
            title('1 night, smoothed nightly ODBA awakenings');
            xticklabels(compose('%1.3f',xticks/3600));
            xlim([1,numel(v)]);
            xlabel('time (h)');
            ylabel('odba');
        end
        for ii = 1:numel(locs)
            allDiffs = [allDiffs,abs(locs - locs(ii))];
        end
    end
    counts = histcounts(allDiffs,linspace(nSmooth,3600*showHours,nBins));
    lns(iPt) = plot(smooth(counts/sum(counts)),'color',[colors(iPt,:),0.2],'linewidth',2);
    hold on;
    allDiffs = [];
end
nTicks = 20;
xticks(linspace(1,numel(counts),nTicks));
xticklabels(compose('%1.1f',linspace(0,showHours,nTicks)));
xtickangle(60);
set(gca,'yscale','log');
ylabel('frac. of awakenings');
xlabel('time (h)');
title({sprintf('distribution of awakenings relative to each other (n=%i)',...
    numel(unique(meta.subjects))),...
    sprintf('ODBA thresh = %1.2f-%1.2f, smoothed/min dist. = %2.0fs',pThreshs(1),pThreshs(end),nSmooth)});
% legend(compose('%1.2f',pThreshs));
legend([lns(1),lns(end)],compose('%1.2f',[pThreshs(1),pThreshs(end)]));
xlim([1 nBins-1]);
ylim([0 10^-1]);


% % % % if strcmp(saveSubject,meta.subjects{iDay}) == 0 ||...
% % % %         (iDay == numel(meta.dates) && strcmp(saveSubject,meta.subjects{end}))
% % % %     saveSubject = meta.subjects{iDay};
% % % %     ff(600,600);
% % % %     histogram(allDiffs,linspace(nSmooth,3600*showHours,nBins));
% % % %     xticklabels(compose('%1.3f',xticks/3600));
% % % %     set(gca,'yscale','log');
% % % %     ylabel('count');
% % % %     xlabel('time (h)');
% % % %     title({sprintf('subject %s',saveSubject),...
% % % %         'distribution of awakenings relative to each other',...
% % % %         sprintf('ODBA thresh = %1.2f, smoothed/min dist. = %2.0fs',pThresh,nSmooth)});
% % % %     allDiffs = [];
% % % % end