% load file
useFiles = {'C1_Jul1_2014__20140623.mat','C0_Jan15_2017__20170116.mat'};
dataPath = '/Users/matt/Documents/Data/KRSP/CompressedAxy';
nBins = 25;
tClock = linspace(0,24,nBins);
colors = lines(numel(useFiles));

close all
ff(800,800,2);
for iFile = 1:numel(useFiles)
    load(fullfile(dataPath,useFiles{iFile}));
    
    dtDays = datetime(year(T.datetime),month(T.datetime),day(T.datetime));
    dtS = secDay(T.datetime);
    allDates = unique(dtDays);
    
    sum_zOdba = NaN(nBins-1,1);
    for iBin = 1:nBins-1
        t_start = tClock(iBin) * 3600;
        t_end = tClock(iBin+1) * 3600;
        dtIds = find(dtS >= t_start & dtS < t_end);
        sum_zOdba(iBin) = sum(T.odba(dtIds));
    end
    polarClock = linspace(0,2*pi,nBins);
    h = polarhistogram('BinEdges',polarClock,'BinCounts',sum_zOdba);
    hold on;
    h.FaceColor = colors(iFile,:);
    h.FaceAlpha = 0.25;
    pax = gca;
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 18;
    pax.Layer = 'top';
    % rlim([0 rlimVal]);
    rticks([]);
    pax.Color = [1 1 1];
    pax.ThetaTick = linspace(0,360,25);
    pax.ThetaTickLabels = compose('%i',0:23);
end
title('Hourly sum(ODBA)');
legend(useFiles,'location','south');