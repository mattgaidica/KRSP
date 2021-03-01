% setup /Users/matt/Documents/MATLAB/KRSP/Figures/figure_rhythmicSleep_noDiff.m

months =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

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
%%
c = colorbar('location','southoutside');
colormap(colors);
c.Limits = [0,1];
c.Ticks = linspace(0,1,12);
c.TickLabels = months;
c.TickDirection = 'out';
title('Asleep');
%%
fs = 14;
text(0,-1,["100%","\downarrow"],'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));
text(0.55,1.15,"Hour of Day",'horizontalalignment','center','verticalalignment','bottom','fontsize',fs,'color',repmat(0.4,[3,1]));


% % monthDoys = linspace(0,366,13);
% % for iMonth = 1:12
% %     monthData = mean(sq_asleep(sq_doys > monthDoys(iMonth) & sq_doys <= monthDoys(iMonth+1),:));
% %     theseAsleep(iMonth,:) = imgaussfilt(monthData,50,'padding','circular');
% % end
% % 
% % close all
% % ff(600,600);
% % polarTime = linspace(-pi,pi,size(theseAsleep,2)) + pi/2;
% % colors = cool(12);
% % for iMonth = 1:12
% %     polarplot(polarTime,round(theseAsleep(iMonth,:)*50)/50,'color',colors(iMonth,:),'linewidth',4);
% %     hold on;
% % end
% % 
% % pax = gca;
% % pax.ThetaTick = linspace(0,360,25);
% % pax.ThetaTickLabel = 0:23;
% % pax.ThetaZeroLocation = 'top';
% % pax.ThetaDir = 'clockwise';
% % pax.FontSize = 18;
% % pax.Layer = 'top';
% % rlim([0 1]);
% % rticks([]);
% % pax.Color = [1 1 1];