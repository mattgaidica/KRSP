if go
    nSweep = 100;
    tws = [3600/2 3600];
    pThreshs = linspace(0.5,0.95,10);
    pThreshs = .75;%[.75,.95];
    startHour = 0;
    endHour = 6;
    fromWakeup = true;
    if fromWakeup
        sweep = fliplr(round(linspace(startHour*3600+1,endHour*3600,nSweep)));
        sw = -1;
    else
        sweep = round(linspace(startHour*3600+1,endHour*3600,nSweep));
        sw = 1;
    end
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
        for iSubject = 1:numel(subjects)
            v_day = sort(compiled_day_odba{iSubject});
            v_night = sort(compiled_night_odba{iSubject});
            for iiN = 1:numel(pThreshs)
                v_dThresh(iSubject,iiN) = v_day(round(numel(v_day)*pThreshs(iiN)));
                v_nThresh(iSubject,iiN) = v_night(round(numel(v_night)*pThreshs(iiN)));
            end
            v_dMean(iSubject) = mean(v_day);
            v_nMean(iSubject) = mean(v_night);
        end
    end
    
    rs = [];
    ps = [];
    for iTw = 1:numel(tws)
        tw = tws(iTw);
        for iPt = 1:numel(pThreshs)
            for iSw = 1:nSweep
                odba_night_counts = [];
                odba_day_counts = [];
                fprintf('iTw %i/%i, iPt %i/%i, iSw %i/%i \n',...
                    iTw,numel(tws),iPt,numel(pThreshs),iSw,nSweep);
                for iDay = 1:numel(meta.dates)
                    iSubject = find(strcmp(meta.subjects{iDay},subjects));
                    
                    odba_day = smooth(meta.odba_day{iDay},nSmooth);
                    odba_day_counts(iDay) = sum(odba_day - v_dMean(iSubject)) / meta.days_length(iDay);
                    if fromWakeup
                        odba_night = smooth(meta.odba_night{iDay}(end-sweep(iSw)-tw:end-sweep(iSw)+1),nSmooth);
                        sw = -1;
                    else
                        odba_night = smooth(meta.odba_night{iDay}(sweep(iSw):sweep(iSw)+tw-1),nSmooth);
                        sw = 1;
                    end
                    [locs, pks] = peakseek(odba_night,minpeakdist,v_nThresh(iSubject,iPt));
                    odba_night_counts(iDay) = numel(locs);
                    % correlated with nightly sum has minimal correlation
%                     odba_night_counts(iDay) = sum(odba_night - v_nMean(iSubject)) / meta.nights_length(iDay);
                end
                [A,I] = rmoutliers(odba_day_counts);
                B = odba_night_counts(~I);
                [B,I] = rmoutliers(B);
                A = A(~I);
                
                [r,pval] = corr(A',B');
                rs(iSw,iPt,iTw) = r;
                ps(iSw,iPt,iTw) = pval;
            end
        end
    end
    go = false;
end

% close all
ff(1200,900);
rows = 5;
cols = size(rs,3);
colors = cool(size(rs,2));
lw = 1;
lns = [];
t = sw*sweep/3600;
for iPt = 1:size(rs,2)
    for iTw = 1:size(rs,3)
        subplot(rows,cols,prc(cols,[1,iTw]));
        rData = squeeze(rs(:,iPt,iTw));
        plot(t,rData,'color',colors(iPt,:),'linewidth',lw);
        hold on;
        ylabel('r');
        title(sprintf('%1.2f hr window',tws(iTw)/3600));
        ylim([-.4 .4]);
        xlim([min(t),max(t)]);
        grid on;
        
        subplot(rows,cols,prc(cols,[2,iTw]));
        pData = squeeze(ps(:,iPt,iTw));
        lns(iPt) = plot(t,pData,'color',colors(iPt,:),'linewidth',lw);
        hold on;
        ylabel('p');
        ylim([0 0.2]);
        xlim([min(t),max(t)]);
        yticks([0,0.01,0.05,0.2]);
        grid on;
    end
end
% legend(lns,compose('%1.2f',pThreshs),'location','northeast');
% legend box off;

rHist = [];
pHist = [];
for iTw = 1:size(rs,3)
    subplot(rows,cols,prc(cols,[3,iTw]));
    rData = [];
    all_m = [];
    for iSw = 1:size(rs,1)
        [m,idx] = max(abs(squeeze(rs(iSw,:,iTw))));
        all_m(iSw) = m;
        rData(iSw) = idx;
    end
    for iPt = 1:size(rs,2)
        plot([min(t),max(t)],[iPt,iPt],'color',[colors(iPt,:),0.5],'linewidth',lw);
        hold on;
    end
    plot(t,rData,'k.','markerSize',15);
    [m_val,m_id] = max(all_m);
    plot(t(m_id),rData(m_id),'ro','markerSize',10);
    text(t(m_id),rData(m_id)+0.2,sprintf('%1.2f',m_val),'color','r',...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
    xlim([min(t),max(t)]);
    ylim([0,size(rs,2)+1]);
    yticks([1:size(rs,2)]);
    yticklabels(compose('%1.2f',pThreshs));
    title('max r');
    grid on;
    ylabel('thresh');
    
    subplot(rows,cols,prc(cols,[4,iTw]));
    pData = [];
    for iSw = 1:size(rs,1)
        [m,idx] = min(squeeze(ps(iSw,:,iTw)));
        all_m(iSw) = m;
        pData(iSw) = idx;
    end
    for iPt = 1:size(rs,2)
        plot([min(t),max(t)],[iPt,iPt],'color',[colors(iPt,:),0.5],'linewidth',lw);
        hold on;
    end
    plot(t,pData,'k.','markerSize',15);
    [m_val,m_id] = min(all_m);
    plot(t(m_id),pData(m_id),'ro','markerSize',10);
    text(t(m_id),pData(m_id)+0.2,sprintf('%1.5f',m_val),'color','r',...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
    xlim([min(t),max(t)]);
    ylim([0,size(rs,2)+1]);
    yticks([1:size(rs,2)]);
    yticklabels(compose('%1.2f',pThreshs));
    title('min p');
    grid on;
    ylabel('thresh');
    if fromWakeup
        xlabel('time from wakeup');
    else
        xlabel('time after going to sleep');
    end
    
    subplot(rows,cols,prc(cols,[5,iTw]));
    rcounts = histcounts(rData,[0:numel(pThreshs)]+0.5);
    pcounts = histcounts(pData,[0:numel(pThreshs)]+0.5);
    bar([rcounts',pcounts']);
    legend('r','p');
    legend box off;
    xticks(1:numel(pThreshs));
    xticklabels(compose('%1.2f',pThreshs));
    xtickangle(90);
    ylabel('count');
    xlabel('thresh');
end