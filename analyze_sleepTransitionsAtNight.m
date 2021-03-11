% setup: /Users/matt/Documents/MATLAB/KRSP/predict_awake.m
% trans_at = [trans_at secDay(Tawake.datetime)'];
% trans_to = [trans_to Tawake.awake'];
% trans_on = [trans_on day(Tawake.datetime,'dayofyear')'];
ssPath = '/Users/matt/Documents/Data/KRSP/SunriseSunset';
files = dir(fullfile(ssPath,'*.txt'));
Tss = readtable(fullfile(ssPath,files(4).name)); % 366 day year (simplify for now)
    
close all
ff(1300,900);
% useDoys = 190:260;
binEdges = linspace(1,86400,100); % rm 0 entries, are those subsequent animal entries?

nHalfWindow = 30;
allDoys = 1:366;
colors = flip(mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png',366));
op = 0.2;
nS = 3;
t = linspace(0,24,numel(binEdges)-1);

for iSun = 1:2
    for iDoy = 1:366
        shiftDoys = circshift(allDoys,-iDoy+1+nHalfWindow);
        useDoys = shiftDoys(1:nHalfWindow*2+1);

        shiftBy = closest(t,secDay(Tss.solar_noon(iDoy))/60/60) + round(numel(t)/2);
        if iSun == 1
            shiftBy = closest(t,secDay(Tss.sunrise(iDoy))/60/60) + round(numel(t)/2);
        end

        subplot(2,2,prc(2,[iSun,1]));
        useIds = trans_to==1 & ismember(trans_on,useDoys);
        counts = histcounts(trans_at(useIds),binEdges) / sum(useIds);
        plot(t,circshift(imgaussfilt(counts,nS,'padding','circular'),-shiftBy),'color',[colors(iDoy,:),op]);
        hold on;
        if iDoy == 1
            title('transition to awake');
            xlim([0 24]);
            ylim([0 0.025]);
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
            shiftBy = closest(t,secDay(Tss.sunset(iDoy))/60/60) + round(numel(t)/2);
        end

        subplot(2,2,prc(2,[iSun,2]));
        useIds = trans_to==0 & ismember(trans_on,useDoys);
        counts = histcounts(trans_at(useIds),binEdges) / sum(useIds);
        plot(t,circshift(imgaussfilt(counts,nS,'padding','circular'),-shiftBy),'color',[colors(iDoy,:),op]);
        hold on;
        if iDoy == 1
            title('transition to sleep');
            xlim([0 24]);
            ylim([0 0.025]);
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

        drawnow;
    end
end

%% CDF plots
sIds = round(linspace(1,366,5));
seasonDoys = circshift(1:366,57); % centers at 171, so light is equal
% seasonDoys = circshift(1:366,57-21); % centers at 192, so temp is equal
useDoys = {seasonDoys(sIds(1):sIds(2)),seasonDoys(sIds(2):sIds(3)),...
    seasonDoys(sIds(3):sIds(4)),seasonDoys(sIds(4):sIds(5))};
colors = flip(mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png',5));

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
            theseDurations = theseTransEnd - theseTransStart;
            % fix transitions around next day (in seconds)
            for ii = find(theseDurations < 0)
                theseDurations(ii) = theseDurations(ii) + 86400;
            end
            sleepDurations = [sleepDurations theseDurations];
        end
    end
    
    % uncomment for PDF
%     [counts,edges] = histcounts(transTimes/3600,linspace(0,8,(8/24)*1440));
%     plot(edges(1:end-1),normalize(smoothdata(counts,'gaussian',20),'range'),'color',colors(iSeason,:),'linewidth',2);
%     hold on; 
    
    % uncomment for CDF
    sleepDurations_season{iSeason} = sleepDurations;
    if ~isempty(sleepDurations)
        [Fout,x,Flo,Fup] = ecdf(sleepDurations/3600);
        lns(iSeason) = plot(x,Fout,'color',colors(iSeason,:),'linewidth',lw1);
        hold on;
        plot(x,Flo,':','color',colors(iSeason,:),'linewidth',lw2);
        plot(x,Fup,':','color',colors(iSeason,:),'linewidth',lw2);
    end
end

xlim([0 showHours]);
set(gca,'fontsize',16);
xlabel('Asleep Length (hrs)');
ylabel('Probability');
grid on;
title('Asleep Periods');
legend(lns,seasonLabels,'location','southeast','AutoUpdate','off');

% do significance
binEdges = linspace(0,86400,24+1);
yText = flip(linspace(0.4,0.8,6));
iCount = 0;
fs = 12;
lw = 2;
pThresh = 0.001;
sigColor = repmat(0.3,[1,3]);
for ii = 1:4
    [~,~,binii] = histcounts(sleepDurations_season{ii},binEdges);
    for jj = ii+1:4
        seasonCompare = sprintf('%s - %s\n',seasonLabels{ii},seasonLabels{jj});
        disp(seasonCompare);
        iCount = iCount + 1;
        text(5,yText(iCount)-.021,seasonCompare,'horizontalalignment','left',...
            'verticalalignment','middle','fontsize',fs,'color',sigColor);
        
        [~,~,binjj] = histcounts(sleepDurations_season{jj},binEdges);
        for iBin = 1:numel(binEdges)-1
            if sum(binii==iBin) > 1 && sum(binjj==iBin) > 1
                x = [sleepDurations_season{ii}(binii==iBin),sleepDurations_season{jj}(binjj==iBin)];
                group = [ones(1,sum(binii==iBin)),zeros(1,sum(binjj==iBin))];
                p = anova1(x,group,'off');
            end
            if p < pThresh && ii ~= jj
                ln = plot([binEdges(iBin)/3600,binEdges(iBin+1)/3600],repmat(yText(iCount),[1,2]),'-',...
                    'linewidth',lw,'color',sigColor);
                uistack(ln,'bottom');
            end
        end
    end
end

% mast plots
useSubplots = [3,4,7,8];
for iSeason = 1:4
    lns = NaN(2,1);
    sleepDurations_mast = {};
    for iMast = 1:2
        if iMast == 1
            useYears = [2015:2018,2020];
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
            theseDurations = theseTransEnd - theseTransStart;
            % fix transitions around next day (in seconds)
            for ii = find(theseDurations < 0)
                theseDurations(ii) = theseDurations(ii) + 86400;
            end
            sleepDurations = [sleepDurations theseDurations];
        end
        
        sleepDurations_mast{iMast} = sleepDurations;
        if ~isempty(sleepDurations)
            [Fout,x,Flo,Fup] = ecdf(sleepDurations/3600);
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
    
    % end of mast loop, do stats
    [~,~,bin1] = histcounts(sleepDurations_mast{1},binEdges);
    [~,~,bin2] = histcounts(sleepDurations_mast{2},binEdges);
    for iBin = 1:numel(binEdges)-1
        if sum(bin1==iBin) > 1 && sum(bin2==iBin) > 1
            x = [sleepDurations_mast{1}(bin1==iBin),sleepDurations_mast{2}(bin2==iBin)];
            group = [ones(1,sum(bin1==iBin)),zeros(1,sum(bin2==iBin))];
            p = anova1(x,group,'off');
            if p < pThresh
                ln = plot([binEdges(iBin)/3600,binEdges(iBin+1)/3600],repmat(mean(yText),[1,2]),'-',...
                    'linewidth',lw,'color',sigColor);
                uistack(ln,'bottom');
            end
        end
    end
end

%% test each year
close all
ff(600,600);
iSeason = 4;
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
        theseDurations = theseTransEnd - theseTransStart;
        % fix transitions around next day (in seconds)
        for ii = find(theseDurations < 0)
            theseDurations(ii) = theseDurations(ii) + 86400;
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
        [Fout,x,Flo,Fup] = ecdf(sleepDurations/3600);
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
title('Asleep Periods, Autumn');
%%

% pd = fitdist(transTimes','Lognormal');
% [m,v] = lognstat(pd.mu,pd.sigma);
% [p,pLo,pUp] = logncdf(transTimes,pd.mu,pd.sigma,pd.ParameterCovariance);
% 
% plot(sort(transTimes),sort(p),'k-');
% hold on;
% plot(sort(transTimes),sort(pLo),'k:');
% plot(sort(transTimes),sort(pUp),'k:');
% 
% [f,x,flo,fup] = ecdf(transTimes);
% plot(x,f,'b');
% plot(x,flo,'b:');
% plot(x,fup,'b:');