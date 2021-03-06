% setup /Users/matt/Documents/MATLAB/KRSP/Figures/figure_rhythmicSleep_noDiff.m
% doing it by months doesn't work, too hard to follow
months =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

%% sleep mean by SEASON
close all
ff(1200,900);
iSubplot = 1;
for iYear = [2014:2017,2019,2020]
    theseAsleep = [];
    for iDoy = 1:366
        useIds = find(ismember(sq_doys,iDoy) & ismember(sq_years,iYear));
        if ~isempty(useIds) & numel(useIds) > 1
            monthData = mean(sq_asleep(useIds,:));
            theseAsleep(iDoy,:) = imgaussfilt(monthData,1440/24,'padding','circular');
        else
            theseAsleep(iDoy,:) = NaN(1,1440);
        end
    end
    
    subplot(2,3,iSubplot);
    iSubplot = iSubplot + 1;
    polarTime = linspace(-pi,pi,size(theseAsleep,2)) + pi/2;
    colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png',366);
    for iDoy = 1:366
        polarplot(polarTime,theseAsleep(iDoy,:),'color',[colors(iDoy,:),0.2],'linewidth',2);
        hold on;
    end

    pax = gca;
    pax.ThetaTick = linspace(0,360,25);
    pax.ThetaTickLabel = 0:23;
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 14;
    pax.Layer = 'top';
    rlim([0 1]);
    rticks([]);
    pax.Color = [1 1 1];

    c = colorbar('location','southoutside');
    colormap(colors);
    c.Limits = [0,1];
    c.Ticks = linspace(0,1,12);
    c.TickLabels = months;
    c.TickDirection = 'out';
    title(sprintf('Asleep Mean %i',iYear));

    fs = 14;
    text(0,-1,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
    text(0.55,1.15,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
    drawnow;
end

%% sleep mean by SEX
ff(1200,900);
iSubplot = 1;
for iSex = 0:1
    for iYear = [2014:2017,2019,2020]
        theseAsleep = [];
        for iDoy = 1:366
            useIds = find(ismember(sq_doys,iDoy) & ismember(sq_years,iYear) & ismember(sq_sex,iSex));
            if ~isempty(useIds) & numel(useIds) > 1
                monthData = mean(sq_asleep(useIds,:));
                theseAsleep(iDoy,:) = imgaussfilt(monthData,1440/24,'padding','circular');
            else
                theseAsleep(iDoy,:) = NaN(1,1440);
            end
        end

        subplot(2,3,iSubplot);
        iSubplot = iSubplot + 1;
        polarTime = linspace(-pi,pi,size(theseAsleep,2)) + pi/2;
        colors = lines(2);
        for iDoy = 1:366
            polarplot(polarTime,theseAsleep(iDoy,:),'color',[colors(iSex+1,:),0.2],'linewidth',2);
            hold on;
        end

        pax = gca;
        pax.ThetaTick = linspace(0,360,25);
        pax.ThetaTickLabel = 0:23;
        pax.ThetaZeroLocation = 'top';
        pax.ThetaDir = 'clockwise';
        pax.FontSize = 14;
        pax.Layer = 'top';
        rlim([0 1]);
        rticks([]);
        pax.Color = [1 1 1];

        c = colorbar('location','southoutside');
        colormap(colors);
        c.Limits = [0,1];
        c.Ticks = linspace(0,1,12);
        c.TickLabels = months;
        c.TickDirection = 'out';
        title(sprintf('Asleep Mean %i',iYear));

        fs = 14;
        text(0,-1,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
        text(0.55,1.15,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
        drawnow;
    end
    iSubplot = 1;
end

% setup /Users/matt/Documents/MATLAB/KRSP/Figures/figure_rhythmicSleep_noDiff.m
% doing it by months doesn't work, too hard to follow
months =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

%% sleep consistency by SEASON
close all
ff(1200,900);
iSubplot = 1;
for iYear = [2014:2017,2019,2020]
    theseAsleep = [];
    for iDoy = 1:366
        useIds = find(ismember(sq_doys,iDoy) & ismember(sq_years,iYear));
        if ~isempty(useIds) & numel(useIds) > 1
            monthData = (1/3)-var(sq_asleep(useIds,:));
            theseAsleep(iDoy,:) = imgaussfilt(monthData,1440/24,'padding','circular');
        else
            theseAsleep(iDoy,:) = NaN(1,1440);
        end
    end
    
    subplot(2,3,iSubplot);
    iSubplot = iSubplot + 1;
    polarTime = linspace(-pi,pi,size(theseAsleep,2)) + pi/2;
    colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png',366);
    for iDoy = 1:366
        polarplot(polarTime,theseAsleep(iDoy,:),'color',[colors(iDoy,:),0.2],'linewidth',2);
        hold on;
    end

    pax = gca;
    pax.ThetaTick = linspace(0,360,25);
    pax.ThetaTickLabel = 0:23;
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 14;
    pax.Layer = 'top';
    rlim([0 1/3]);
    rticks([]);
    pax.Color = [1 1 1];

    c = colorbar('location','southoutside');
    colormap(colors);
    c.Limits = [0,1];
    c.Ticks = linspace(0,1,12);
    c.TickLabels = months;
    c.TickDirection = 'out';
    title(sprintf('Asleep Consistency %i',iYear));

    fs = 14;
    text(0,-1/3,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
    text(0.55,.38,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
    drawnow;
end

%% sleep consistency by SEX
ff(1200,900);
iSubplot = 1;
for iSex = 0:1
    for iYear = [2014:2017,2019,2020]
        theseAsleep = [];
        for iDoy = 1:366
            useIds = find(ismember(sq_doys,iDoy) & ismember(sq_years,iYear) & ismember(sq_sex,iSex));
            if ~isempty(useIds) & numel(useIds) > 1
                monthData = (1/3)-var(sq_asleep(useIds,:));
                theseAsleep(iDoy,:) = imgaussfilt(monthData,1440/24,'padding','circular');
            else
                theseAsleep(iDoy,:) = NaN(1,1440);
            end
        end

        subplot(2,3,iSubplot);
        iSubplot = iSubplot + 1;
        polarTime = linspace(-pi,pi,size(theseAsleep,2)) + pi/2;
        colors = lines(2);
        for iDoy = 1:366
            polarplot(polarTime,theseAsleep(iDoy,:),'color',[colors(iSex+1,:),0.2],'linewidth',2);
            hold on;
        end

        pax = gca;
        pax.ThetaTick = linspace(0,360,25);
        pax.ThetaTickLabel = 0:23;
        pax.ThetaZeroLocation = 'top';
        pax.ThetaDir = 'clockwise';
        pax.FontSize = 14;
        pax.Layer = 'top';
        rlim([0 1/3]);
        rticks([]);
        pax.Color = [1 1 1];

        c = colorbar('location','southoutside');
        colormap(colors);
        c.Limits = [0,1];
        c.Ticks = linspace(0,1,12);
        c.TickLabels = months;
        c.TickDirection = 'out';
        title(sprintf('Asleep Consistency %i',iYear));

        fs = 14;
        text(0,-1/3,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
        text(0.55,.38,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
        drawnow;
    end
    iSubplot = 1;
end