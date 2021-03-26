% setup: /Users/matt/Documents/MATLAB/KRSP/predict_awake.m
% trans_at = [trans_at secDay(Tawake.datetime)'];
% trans_to = [trans_to Tawake.awake'];
% trans_on = [trans_on day(Tawake.datetime,'dayofyear')'];

%% this is all the transitions (x-y plot) for all seasons
% sunrise/set centered on top, solar noon on bottom
Tss = makeTss(2014:2020);
months =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

close all
ff(1300,900);
binEdges = linspace(1,86400,100); % rm 0 entries, are those subsequent animal entries?

nHalfWindow = 30;
allDoys = 1:366;
colors = seasonColors(1:366); % this will use all doys because it's windowed
op = 0.2;
nS = 3;
t = linspace(0,24,numel(binEdges)-1);

useylim = 0.03;
for iSun = 1:2
    for iDoy = 1:366
        shiftDoys = circshift(allDoys,-iDoy+1+nHalfWindow);
        useDoys = shiftDoys(1:nHalfWindow*2+1);
        shiftBy = closest(t,mean(secDay(Tss.noon(Tss.doy == iDoy)),1)/60/60) + round(numel(t)/2);
        if iSun == 1
            shiftBy = closest(t,mean(secDay(Tss.sunrise(Tss.doy == iDoy)),1)/60/60) + round(numel(t)/2);
        end

        subplot(2,2,prc(2,[iSun,1]));
        useIds = trans_to==1 & ismember(trans_on,useDoys);
        counts = histcounts(trans_at(useIds),binEdges) / sum(useIds);
        plot(t,circshift(imgaussfilt(counts,nS,'padding','circular'),-shiftBy),'color',[colors(iDoy,:),op]);
        hold on;
        if iDoy == 1
            title('transition to awake');
            xlim([0 24]);
            ylim([0 useylim]);
            ylabel('probability')
            if iSun == 1
                xlabel('hours relative to sunrise');
            else
                xlabel('hours relative to solar noon');
            end
            xticks(0:6:24);
            xticklabels({'±12','-6','0','+6','±12'});
            set(gca,'fontsize',16)
            grid on;
            c = colorbar('location','southoutside');
            colormap(colors);
            c.Limits = [0,1];
            c.Ticks = linspace(0,1,12);
            c.TickLabels = months;
            c.TickDirection = 'out';
            c.FontSize = 9;
            plot([12 12],ylim,'k:');
            plot([36 36],ylim,'k:');
        end

        if iSun == 1
            shiftBy = closest(t,mean(secDay(Tss.sunset(Tss.doy == iDoy)),1)/60/60) + round(numel(t)/2);
        end

        subplot(2,2,prc(2,[iSun,2]));
        useIds = trans_to==0 & ismember(trans_on,useDoys);
        counts = histcounts(trans_at(useIds),binEdges) / sum(useIds);
        plot(t,circshift(imgaussfilt(counts,nS,'padding','circular'),-shiftBy),'color',[colors(iDoy,:),op]);
        hold on;
        if iDoy == 1
            title('transition to sleep');
            xlim([0 24]);
            ylim([0 useylim]);
            xticks(0:6:24);
            xticklabels({'±12','-6','0','+6','±12'});
            ylabel('probability');
            if iSun == 1
                xlabel('hours relative to sunset');
            else
                xlabel('hours relative to solar noon');
            end
            set(gca,'fontsize',16)
            grid on;
            c = colorbar('location','southoutside');
            colormap(colors);
            c.Limits = [0,1];
            c.Ticks = linspace(0,1,12);
            c.TickLabels = months;
            c.TickDirection = 'out';
            c.FontSize = 9;
            plot([12 12],ylim,'k:');
            plot([36 36],ylim,'k:');
        end
    end
    drawnow;
end

%% CDF plots
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);

% main plot
close all
ff(1000,400);
lw1 = 3;
lw2 = 1;
rows = 2;
cols = 4;
subplot(rows,cols,[1 2 5 6]);
lns = [];
showHours = 8;
seasonLabels = {'Winter','Spring','Summer','Autumn'};
sleepDurations_season = {};
for iSeason = 1:4
    uniqueSqs = unique(trans_is);
    sleepDurations = [];
    for iSq = uniqueSqs
        theseSleepTrans = find(trans_is==iSq & trans_to==0 & ismember(trans_on,useDoys{iSeason}));
        if ~isempty(theseSleepTrans)
            
            % test light/dark transitions
% % % %             theseDoys = trans_on(theseSleepTrans);
% % % %             meanSunrise = mean(secDay(Tss.sunrise(theseDoys)));
% % % %             meanSunset = mean(secDay(Tss.sunset(theseDoys)));
% % % %             temp = [];
% % % %             for ii = theseSleepTrans
% % % %                 if trans_at(ii) > meanSunset || trans_at(ii) < meanSunrise
% % % %                     temp = [temp ii];
% % % %                 end
% % % %             end
% % % %             if isempty(temp)
% % % %                 continue;
% % % %             end
% % % %             theseSleepTrans = temp;
            
            allEntries = find(trans_is==iSq);
            % must end with awake transition to get time
            if trans_to(allEntries(end)) == 0
                theseSleepTrans = theseSleepTrans(1:end-1);
            end

            theseTransStart = trans_at(theseSleepTrans); % asleep
            theseTransEnd = trans_at(theseSleepTrans+1); % awake
            % make minutes, that's the highest resolution from
            % predict_awake anyways
            theseDurations = (theseTransEnd - theseTransStart) / 60;
            % fix transitions around next day (in seconds)
            for ii = find(theseDurations < 0)
                theseDurations(ii) = theseDurations(ii) + 1440;
            end
            sleepDurations = [sleepDurations theseDurations];
        end
    end
    
    % uncomment for PDF
%     [counts,edges] = histcounts(sleepDurations/60,linspace(0,8,(8/24)*1440));
%     plot(edges(1:end-1),normalize(smoothdata(counts,'gaussian',5),'range'),'color',colors(iSeason,:),'linewidth',2);
%     hold on; 
    
    % uncomment for CDF
    sleepDurations_season{iSeason} = sleepDurations;
    if ~isempty(sleepDurations)
        [Fout,x,Flo,Fup] = ecdf(sleepDurations/60); % hours
        lns(iSeason) = plot(x,Fout,'color',colors(iSeason,:),'linewidth',lw1);
        hold on;
        plot(x,Flo,':','color',colors(iSeason,:),'linewidth',lw2);
        plot(x,Fup,':','color',colors(iSeason,:),'linewidth',lw2);
    end
end

% do stats
ylab = flip(linspace(0.4,0.8,6));
statCount = 0;
fs = 12;
statColor = repmat(0.5,[1,4]);
for ii = 1:4
    for jj = ii+1:4
        % if p < 0.05 dist. are different
        [h,p] = kstest2(sleepDurations_season{ii},sleepDurations_season{jj});
        sigStr = '';
        if p < 0.05 && p >= 0.01
            sigStr = '*';
        elseif p < 0.01 && p >= 0.001
            sigStr = '**';
        elseif p < 0.001
            sigStr = '***';
        end
        statCount = statCount + 1;
        text(4,ylab(statCount),sprintf('%s-%s %sp = %1.2e',seasonLabels{ii}(1:3),...
            seasonLabels{jj}(1:3),sigStr,p),'fontsize',fs,'color',statColor);
    end
end


xlim([0 showHours]);
set(gca,'fontsize',16);
xlabel('Asleep Length (hrs)');
ylabel('Probability');
grid on;
title('Asleep Periods');
legend(lns,seasonLabels,'location','southeast','AutoUpdate','off');

% mast plots
useSubplots = [3,4,7,8];
for iSeason = 1:4
    lns = NaN(2,1);
    sleepDurations_mast = {};
    for iMast = 1:2
        if iMast == 1
            useYears = [2015:2020];
            op = 1;
        else
            useYears = [2014,2019];
            op = 0.4;
        end
        subplot(rows,cols,useSubplots(iSeason));
        uniqueSqs = unique(trans_is);
        sleepDurations = [];
        for iSq = uniqueSqs
            theseSleepTrans = find(trans_is==iSq & trans_to==0 &...
                ismember(trans_on,useDoys{iSeason}) & ismember(trans_yr,useYears));
            
            % test light/dark transitions
% % % %             theseDoys = trans_on(theseSleepTrans);
% % % %             meanSunrise = mean(secDay(Tss.sunrise(theseDoys)));
% % % %             meanSunset = mean(secDay(Tss.sunset(theseDoys)));
% % % %             temp = [];
% % % %             for ii = theseSleepTrans
% % % %                 if trans_at(ii) > meanSunset || trans_at(ii) < meanSunrise
% % % %                     temp = [temp ii];
% % % %                 end
% % % %             end
% % % %             if isempty(temp)
% % % %                 continue;
% % % %             end
% % % %             theseSleepTrans = temp;

            allEntries = find(trans_is==iSq);
            % must end with awake transition to get time
            if trans_to(allEntries(end)) == 0
                theseSleepTrans = theseSleepTrans(1:end-1);
            end

            theseTransStart = trans_at(theseSleepTrans); % asleep
            theseTransEnd = trans_at(theseSleepTrans+1); % awake
            theseDurations = (theseTransEnd - theseTransStart) / 60;
            % fix transitions around next day (in seconds)
            for ii = find(theseDurations < 0)
                theseDurations(ii) = theseDurations(ii) + 1440;
            end
            sleepDurations = [sleepDurations theseDurations];
        end
        
        sleepDurations_mast{iMast} = sleepDurations;
        if ~isempty(sleepDurations)
            [Fout,x,Flo,Fup] = ecdf(sleepDurations/60);
            lns(iMast) = plot(x,Fout,'color',[colors(iSeason,:),op],'linewidth',lw1);
            hold on;
            plot(x,Flo,':','color',[colors(iSeason,:),op],'linewidth',lw2);
            plot(x,Fup,':','color',[colors(iSeason,:),op],'linewidth',lw2);
        end
        xlim([0 showHours]);
        set(gca,'fontsize',12);
        xlabel('Asleep Length (hrs)');
        grid on;
        title(seasonLabels{iSeason});
        if ~isnan(lns(1)) && ~isnan(lns(2))
            legend(lns,{'Non-mast','Mast'},'location','southeast','AutoUpdate','off');
        elseif ~isnan(lns(1))
            legend(lns(1),{'Non-mast'},'location','southeast','AutoUpdate','off');
        else
            legend(lns(2),{'Mast'},'location','southeast','AutoUpdate','off');
        end
    end
    
    if ~isempty(sleepDurations_mast{1}) && ~isempty(sleepDurations_mast{2})
        [h,p] = kstest2(sleepDurations_mast{1},sleepDurations_mast{2});
        sigStr = '';
        if p < 0.05 && p >= 0.01
            sigStr = '*';
        elseif p < 0.01 && p >= 0.001
            sigStr = '**';
        elseif p < 0.001
            sigStr = '***';
        end
            
        text(2,mean(ylab),sprintf('%sp = %1.2e',sigStr,p),'fontsize',fs,'color',statColor);
    else
        text(2,mean(ylab),'no data','fontsize',fs,'color',statColor);
    end
end

%% test each year
seasonLabels = {'Winter','Spring','Summer','Autumn'};
lw1 = 3;
lw2 = 1;
close all
ff(1200,900);
% iSeason = 4;
for iSeason = 1:4
    subplot(2,2,iSeason);
    lns = [];
    legendLabels = {};
    for iYear = 2014:2020
        if ismember(iYear,[2014,2019])
            useColor = 'r';
        else
            useColor = 'k';
        end
        uniqueSqs = unique(trans_is);
        sleepDurations = [];
        for iSq = uniqueSqs
            theseSleepTrans = find(trans_is==iSq & trans_to==0 &...
                ismember(trans_on,useDoys{iSeason}) & trans_yr == iYear);
            nSq = numel(unique(trans_is(ismember(trans_on,useDoys{iSeason}) & trans_yr == iYear)));

            allEntries = find(trans_is==iSq);
            % must end with awake transition to get time
            if trans_to(allEntries(end)) == 0
                theseSleepTrans = theseSleepTrans(1:end-1);
            end

            theseTransStart = trans_at(theseSleepTrans); % asleep
            theseTransEnd = trans_at(theseSleepTrans+1); % awake
            theseDurations = (theseTransEnd - theseTransStart) / 60;
            % fix transitions around next day (in seconds)
            for ii = find(theseDurations < 0)
                theseDurations(ii) = theseDurations(ii) + 1440;
            end
            sleepDurations = [sleepDurations theseDurations];
        end
        % uncomment for PDF
    %     [counts,edges] = histcounts(transTimes/3600,linspace(0,8,(8/24)*1440));
    %     plot(edges(1:end-1),normalize(smoothdata(counts,'gaussian',20),'range'),'color',colors(iSeason,:),'linewidth',2);
    %     hold on; 

        % uncomment for CDF
        if ~isempty(sleepDurations)
            legendLabels{numel(legendLabels)+1} = sprintf('%i, n = %i',iYear,nSq);
            [Fout,x,Flo,Fup] = ecdf(sleepDurations/60);
            lns(numel(lns)+1) = plot(x,Fout,'color',useColor,'linewidth',lw1);
            hold on;
            plot(x,Flo,':','color',useColor,'linewidth',lw2);
            plot(x,Fup,':','color',useColor,'linewidth',lw2);
        end
    end
    legend(lns,legendLabels,'location','southeast');
    xlim([0 8]);
    set(gca,'fontsize',16);
    xlabel('Asleep Length (hrs)');
    ylabel('%');
    grid on;
    title(sprintf('Asleep Periods, %s',seasonLabels{iSeason}));
end