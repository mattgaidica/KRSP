% setup /Users/matt/Documents/MATLAB/KRSP/Figures/figure_rhythmicSleep_noDiff.m
% doing it by months doesn't work, too hard to follow
months =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

%% sleep mean
theseAsleep = [];
for iDoy = 1:366
    if sum(sq_doys == iDoy) > 1
        monthData = mean(sq_asleep(sq_doys == iDoy,:));
        theseAsleep(iDoy,:) = imgaussfilt(monthData,1440/24,'padding','circular');
    else
        theseAsleep(iDoy,:) = NaN(1,1440);
    end
end

close all
ff(500,500);
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
title('Asleep Mean');

fs = 14;
text(0,-1,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
text(0.55,1.15,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));

%% sleep consistency
theseAsleep = [];
for iDoy = 1:366
    if sum(sq_doys == iDoy) > 4
        monthData = (1/3)-var(sq_asleep(sq_doys == iDoy,:));
        theseAsleep(iDoy,:) = imgaussfilt(monthData,1440/24,'padding','circular');
    else
        theseAsleep(iDoy,:) = NaN(1,1440);
    end
end

% close all
ff(500,500);
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
title('Asleep Consistency');

fs = 14;
text(0,-1/3,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
text(0.55,.38,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));