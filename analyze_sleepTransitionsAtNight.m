% setup: /Users/matt/Documents/MATLAB/KRSP/predict_awake.m
% trans_at = [trans_at secDay(Tawake.datetime)'];
% trans_to = [trans_to Tawake.awake'];
% trans_on = [trans_on day(Tawake.datetime,'dayofyear')'];
ssPath = '/Users/matt/Documents/Data/KRSP/SunriseSunset';
files = dir(fullfile(ssPath,'*.txt'));
T_ss = readtable(fullfile(ssPath,files(4).name)); % 366 day year (simplify for now)
    
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

        shiftBy = closest(t,secDay(T_ss.solar_noon(iDoy))/60/60) + round(numel(t)/2);
        if iSun == 1
            shiftBy = closest(t,secDay(T_ss.sunrise(iDoy))/60/60) + round(numel(t)/2);
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
            shiftBy = closest(t,secDay(T_ss.sunset(iDoy))/60/60) + round(numel(t)/2);
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
seasonLabels = {'Winter','Spring','Summer','Autumn'};
for iSeason = 1:4
    uniqueSqs = unique(trans_is);
    transTimes = [];
    for iSq = uniqueSqs
        theseSleepTrans = find(trans_is==iSq & trans_to==0 & ismember(trans_on,useDoys{iSeason}));

        allEntries = find(trans_is==iSq);
        % must end with awake transition to get time
        if trans_to(allEntries(end)) == 0
            theseSleepTrans = theseSleepTrans(1:end-1);
        end

        theseTransStart = trans_at(theseSleepTrans); % asleep
        theseTransEnd = trans_at(theseSleepTrans+1); % awake
        theseTransTimes = theseTransEnd - theseTransStart;
        % fix transitions around next day (in seconds)
        for ii = find(theseTransTimes < 0)
            theseTransTimes(ii) = theseTransTimes(ii) + 86400;
        end
        transTimes = [transTimes theseTransTimes];
    end
    % uncomment for PDF
%     [counts,edges] = histcounts(transTimes/3600,linspace(0,8,(8/24)*1440));
%     plot(edges(1:end-1),normalize(smoothdata(counts,'gaussian',20),'range'),'color',colors(iSeason,:),'linewidth',2);
%     hold on; 
    
    % uncomment for CDF
    [Fout,x,Flo,Fup] = ecdf(transTimes/3600);
    lns(iSeason) = plot(x,Fout,'color',colors(iSeason,:),'linewidth',lw1);
    hold on;
    plot(x,Flo,':','color',colors(iSeason,:),'linewidth',lw2);
    plot(x,Fup,':','color',colors(iSeason,:),'linewidth',lw2);
end
legend(lns,seasonLabels,'location','southeast');
xlim([0 8]);
set(gca,'fontsize',16);
xlabel('Asleep Length (hrs)');
ylabel('%');
grid on;
title('Asleep Periods');

% mast plots
useSubplots = [3,4,7,8];
for iSeason = 1:4
    lns = NaN(2,1);
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
        transTimes = [];
        for iSq = uniqueSqs
            theseSleepTrans = find(trans_is==iSq & trans_to==0 &...
                ismember(trans_on,useDoys{iSeason}) & ismember(trans_yr,useYears));

            allEntries = find(trans_is==iSq);
            % must end with awake transition to get time
            if trans_to(allEntries(end)) == 0
                theseSleepTrans = theseSleepTrans(1:end-1);
            end

            theseTransStart = trans_at(theseSleepTrans); % asleep
            theseTransEnd = trans_at(theseSleepTrans+1); % awake
            theseTransTimes = theseTransEnd - theseTransStart;
            % fix transitions around next day (in seconds)
            for ii = find(theseTransTimes < 0)
                theseTransTimes(ii) = theseTransTimes(ii) + 86400;
            end
            transTimes = [transTimes theseTransTimes];
        end
        
        if ~isempty(transTimes)
            [Fout,x,Flo,Fup] = ecdf(transTimes/3600);
            lns(iMast) = plot(x,Fout,'color',[colors(iSeason,:),op],'linewidth',lw1);
            hold on;
            plot(x,Flo,':','color',[colors(iSeason,:),op],'linewidth',lw2);
            plot(x,Fup,':','color',[colors(iSeason,:),op],'linewidth',lw2);
        end
        xlim([0 8]);
        set(gca,'fontsize',12);
        xlabel('Asleep Length (hrs)');
        ylabel('%');
        grid on;
        title(seasonLabels{iSeason});
        if ~isnan(lns(1)) && ~isnan(lns(2))
            legend(lns,{'Non-mast','Mast'},'location','southeast');
        elseif ~isnan(lns(1))
            legend(lns(1),{'Non-mast'},'location','southeast');
        else
            legend(lns(2),{'Mast'},'location','southeast');
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
    transTimes = [];
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
        theseTransTimes = theseTransEnd - theseTransStart;
        % fix transitions around next day (in seconds)
        for ii = find(theseTransTimes < 0)
            theseTransTimes(ii) = theseTransTimes(ii) + 86400;
        end
        transTimes = [transTimes theseTransTimes];
    end
    % uncomment for PDF
%     [counts,edges] = histcounts(transTimes/3600,linspace(0,8,(8/24)*1440));
%     plot(edges(1:end-1),normalize(smoothdata(counts,'gaussian',20),'range'),'color',colors(iSeason,:),'linewidth',2);
%     hold on; 
    
    % uncomment for CDF
    if ~isempty(transTimes)
        legendLabels{numel(legendLabels)+1} = sprintf('%i, n = %i',iYear,nSq);
        [Fout,x,Flo,Fup] = ecdf(transTimes/3600);
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