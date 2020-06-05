if go
    nSweep = 200;
    tw = 3600/2;
    pThreshs = .75;
    startHour = 0;
    endHour = 6;
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
        v_nMean = [];
        v_dStd = [];
        v_nStd = [];
        for iSubject = 1:numel(subjects)
            v_day = sort(compiled_day_odba{iSubject});
            v_night = sort(compiled_night_odba{iSubject});
            v_dThresh(iSubject) = v_day(round(numel(v_day)*pThreshs(1)));
            v_nThresh(iSubject) = v_night(round(numel(v_night)*pThreshs(1)));
            v_dMean(iSubject) = mean(v_day);
            v_nMean(iSubject) = mean(v_night);
            v_dStd(iSubject) = std(v_day);
            v_nStd(iSubject) = std(v_night);
        end
    end
    
    rs = [];
    ps = [];
    rs_o = [];
    ps_o = [];
    for iWake = 1:2
        if iWake == 1
            sweep = round(linspace(startHour*3600+1,endHour*3600,nSweep));
        else
            sweep = fliplr(round(linspace(startHour*3600+1,endHour*3600,nSweep)));
        end
        for iSw = 1:nSweep
            odba_night_counts = [];
            odba_day_counts = [];
            fprintf('iSw %i/%i \n',iSw,nSweep);
            for iDay = 1:numel(meta.dates)
                iSubject = find(strcmp(meta.subjects{iDay},subjects));
                
                odba_day = smooth(meta.odba_day{iDay},nSmooth);
                odba_day_counts(iDay) = mean((odba_day - v_dMean(iSubject)) / v_dStd(iSubject));% / meta.days_length(iDay);
                if iWake == 1
                    odba_night = smooth(meta.odba_night{iDay}(sweep(iSw):sweep(iSw)+tw-1),nSmooth);
                else
                    odba_night = smooth(meta.odba_night{iDay}(end-sweep(iSw)-tw:end-sweep(iSw)+1),nSmooth);
                end
                [locs, pks] = peakseek(odba_night,minpeakdist,v_nThresh(iSubject));
                odba_night_counts(iDay) = numel(locs);
                % correlated with nightly sum has minimal correlation
%                 odba_night_counts(iDay) = sum(odba_night - v_nMean(iSubject));% / meta.nights_length(iDay);
            end
            % w/ outliers
            [r,p] = corr(odba_day_counts',odba_night_counts');
            rs_o(iWake,iSw) = r;
            ps_o(iWake,iSw) = p;
            % rmoutliers
            [A,I] = rmoutliers(odba_day_counts);
            B = odba_night_counts(~I);
            [B,I] = rmoutliers(B);
            A = A(~I);
            [r,p] = corr(A',B');
            rs(iWake,iSw) = r;
            ps(iWake,iSw) = p;
        end
    end
    go = false;
end

close all;
ylimVals = [-0.35 0.35];
sigThresh = 0.05;
yArrow = ylimVals(2) - 0.05;
lw = 2;
lns = [];
rows = 2;
cols = 1;
t = sweep/3600;
xText = {'time after going to sleep (hrs)','time from wakeup (hrs)'};
fontSize = 14;

h = ff(1200,800);
for iWake = 1:2
    subplot(rows,cols,iWake);
    if iWake == 2
        t = -t;
    end
    lns(1) = plot(t,rs(iWake,:),'linewidth',lw,'color','k');
    hold on;
    lns(2) = plot(t,rs_o(iWake,:),':','linewidth',1,'color','r');
    ylim([ylimVals(1) ylimVals(2)+0.05]);
    yticks(sort([ylimVals,0,ylimVals(2)+0.05]));
    yticklabels({ylimVals(1),0,ylimVals(2),sprintf('p < %1.2f',sigThresh)});
    sig_rs = find(ps(iWake,:) < sigThresh);
    if ~isempty(sig_rs)
        plot(t(sig_rs),ylimVals(2)+0.035,'.','markersize',15,'color','k');
    end
    sig_rs_o = find(ps_o(iWake,:) < sigThresh);
    if ~isempty(sig_rs_o)
        plot(t(sig_rs_o),ylimVals(2)+0.015,'.','markersize',15,'color','r');
    end
    legend(lns,{'rm outliers','w/ outliers'},'location','southeast');
    legend box off;
    grid on;
    box off;
    xlabel(xText{iWake});
    ylabel('r (night awakenings - day ODBA)');
    if iWake == 1
        title({'Do night awakenings affect day odba?',...
            sprintf('awakening > %1.2f subject mean, sliding window %1.2f hrs',pThreshs,tw/3600)});
        text(min(t),yArrow,'\leftarrow squirrel goes to sleep','fontsize',fontSize);
    else
        text(max(t),yArrow,'squirrel wakes up \rightarrow','horizontalalignment','right','fontsize',fontSize);
    end
    set(gca,'fontsize',fontSize);
end