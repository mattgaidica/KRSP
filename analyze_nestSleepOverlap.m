filespath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
if do
    files = dir(fullfile(filespath,'*.mat'));
    overlapStats = [];
    iCount = 0;
    mean_doys = [];
    for iFile = 1:numel(files)
        load(fullfile(filespath,files(iFile).name)); % T, Tstat
        if ~isValidT(T,true)
            disp(['Skipping ',files(iFile).name]);
            continue;
        end
        iCount = iCount + 1;
        T = detect_sleepWake(T,2);
        T.nest_bin = strcmp(T.nest,'Nest');
        overlapStats(iCount,1) = sum(T.nest_bin & T.awake) / size(T,1); % in-awake
        overlapStats(iCount,2) = sum(T.nest_bin & ~T.awake) / size(T,1); % in-asleep
        overlapStats(iCount,3) = sum(~T.nest_bin & T.awake) / size(T,1); % out-awake
        overlapStats(iCount,4) = sum(~T.nest_bin & ~T.awake) / size(T,1); % out-asleep
        mean_doys(iCount) = mean(unique(day(T.datetime,'dayofyear')));
    end
    do = false;
end

if false
    bee_stats = overlapStats(:);
    bee_cats = repmat(1:4,[size(bee_stats,1)/4,1]);
    bee_cats = bee_cats(:);
    catLabels = {'In & Awake','In & Asleep','Out & Awake','Out & Asleep'};
    
    close all
    ff(500,500);
    beeswarm(bee_cats,bee_stats,'corral_style','omit','overlay_style','sd');
    axis tight
    xticks([1:4]);
    xticklabels(catLabels);
    ylabel('frac. of time');
end

close all
% statsMult = [1 -1;-1 1;1 1;-1 -1];
statsMult = [0 -1;0 1;1 0;-1 0];
ff(1800,700);
op = 0.075;
useIds = [2,3,1,4,2];
colors = lines(4);
for iPlot = 1:2
    subplot(1,2,iPlot);
    lns = [];
    for ii = 1:size(overlapStats,1)
        theseStats = overlapStats(ii,:);
        for jj = 1:4
            xs = [theseStats(useIds(jj))*statsMult(useIds(jj),1),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),1)];
            ys = [theseStats(useIds(jj))*statsMult(useIds(jj),2),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),2)];
            lns(1) = plot(xs,ys,'color',[0,0,0,op]);
            hold on;
        end
    end
    if iPlot == 1
        theseStats = median(overlapStats);
        for jj = 1:4
            xs = [theseStats(useIds(jj))*statsMult(useIds(jj),1),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),1)];
            ys = [theseStats(useIds(jj))*statsMult(useIds(jj),2),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),2)];
            lns(2) = plot(xs,ys,'color',[1,0,0,.8],'linewidth',2);
            hold on;
        end
        legend(lns,{'Individual','Median'});
        title('All Data');
    else % seasons
        seasonIds = round(linspace(1,366,5));
        for iSeason = 1:4
            sIds = find(mean_doys >= seasonIds(iSeason) & mean_doys < seasonIds(iSeason+1));
            theseStats = median(overlapStats(sIds,:));
            for jj = 1:4
                xs = [theseStats(useIds(jj))*statsMult(useIds(jj),1),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),1)];
                ys = [theseStats(useIds(jj))*statsMult(useIds(jj),2),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),2)];
                lns(iSeason+1) = plot(xs,ys,'color',[colors(iSeason,:),1],'linewidth',2);
                hold on;
            end
        end
        legend(lns,{'Individual','Winter','Spring','Summer','Fall'});
        title('Seasonal Data');
    end
    xlim([-0.25 1]);
    ylim(xlim);
    xticks([-1 0 1]);
    yticks(xticks);
    xticklabels(abs(xticks));
    yticklabels(xticklabels);
    offset = 0.02;
    text(0,max(ylim)-offset,'in-asleep','horizontalalignment','center','fontsize',12);
    text(max(xlim)-offset,0,'out-awake','horizontalalignment','right','fontsize',12);
    text(0,min(ylim)+offset,'in-awake','horizontalalignment','center','fontsize',12);
    text(min(xlim)+offset,0,'out-asleep','horizontalalignment','left','fontsize',12);
    grid on;
    xlabel('fraction of day');
    ylabel('fraction of day');
    set(gca,'fontsize',14);
end