% setup: /Users/matt/Documents/MATLAB/KRSP/predict_awake.m
% trans_at = [trans_at secDay(Tawake.datetime)'];
% trans_to = [trans_to Tawake.awake'];
% trans_on = [trans_on day(Tawake.datetime,'dayofyear')'];
ssPath = '/Users/matt/Documents/Data/KRSP/SunriseSunset';
files = dir(fullfile(ssPath,'*.txt'));
T_ss = readtable(fullfile(ssPath,files(4).name)); % 366 day year (simplify for now)
    
close all
ff(1300,600);
% useDoys = 190:260;
binEdges = linspace(1,86400,100); % rm 0 entries

nHalfWindow = 30;
allDoys = 1:366;
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png',366);
op = 0.2;
nS = 3;
t = linspace(0,24,numel(binEdges)-1);

for iDoy = 1:366
    shiftDoys = circshift(allDoys,-iDoy+1+nHalfWindow);
    useDoys = shiftDoys(1:nHalfWindow*2+1);
    
    shiftBy = closest(t,secDay(T_ss.sunrise(iDoy))/60/60) + round(numel(t)/2);
    
    subplot(121);
    useIds = trans_to==1 & ismember(trans_on,useDoys);
    counts = histcounts(trans_at(useIds),binEdges) / sum(useIds);
    plot(t,circshift(imgaussfilt(counts,nS,'padding','circular'),-shiftBy),'color',[colors(iDoy,:),op]);
    hold on;
    if iDoy == 1
        title('transition to awake');
        xlim([0 24]);
        ylim([0 0.025]);
        ylabel('probability')
        xlabel('hours relative to sunrise');
        xticks(0:6:24);
        xticklabels({'-12','-6','0','+6','+12'});
        set(gca,'fontsize',16)
        grid on;
        c = colorbar('location','southoutside');
        colormap(colors);
        c.Limits = [0,1];
        c.Ticks = linspace(0,1,12);
        c.TickLabels = months;
        c.TickDirection = 'out';
        c.FontSize = 9;
    end
    
    shiftBy = closest(t,secDay(T_ss.sunset(iDoy))/60/60) + round(numel(t)/2);
    
    subplot(122);
    useIds = trans_to==0 & ismember(trans_on,useDoys);
    counts = histcounts(trans_at(useIds),binEdges) / sum(useIds);
    plot(t,circshift(imgaussfilt(counts,nS,'padding','circular'),-shiftBy),'color',[colors(iDoy,:),op]);
    hold on;
    if iDoy == 1
        title('transition to sleep');
        xlim([0 24]);
        ylim([0 0.025]);
        xticks(0:6:24);
        xticklabels({'-12','-6','0','+6','+12'});
        ylabel('probability');
        xlabel('hours relative to sunset');
        set(gca,'fontsize',16)
        grid on;
        c = colorbar('location','southoutside');
        colormap(colors);
        c.Limits = [0,1];
        c.Ticks = linspace(0,1,12);
        c.TickLabels = months;
        c.TickDirection = 'out';
        c.FontSize = 9;
    end
    
    drawnow;
end