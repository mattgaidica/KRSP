if ~exist('sunYear')
    [sunYear,sunYearAvg,sunMonth,sunHeader,monthNames] = getSunData(1:12);
    [squirrelArr,dailyMean,dailyStd,squirrelHeader] = getSquirrelData();
    yearlyWeather = getWeatherData();
end

colors = lines(8);
sunColor = colors(3,:);
accelColor = colors(5,:); %colors(5,:);
seasonColor = colors(2,:);
mutedYl = [0 0 0];%[0.55 0.43 0.13];
linewidth = 2;
rows = 1;
cols = 2;
close all;
ff(600,600);
rlimVal = 1.55;

% maybe later, use circular color scheme?
% % % % % season outlines
% % % % edges = linspace(0,2*pi,13);
% % % % counts = [1 0 1 1 0 0 1 1 0 0 0 1] * rlimVal;
% % % % h = polarhistogram('BinEdges',edges,'BinCounts',counts,...
% % % %     'FaceColor',seasonColor,'LineWidth',0.5,'FaceAlpha',0.5);
% % % % h.EdgeColor = 'none'; %accelColor;
% % % % hold on;
% % % % % overlay white to make season sectors
% % % % edges = linspace(0,2*pi,13);
% % % % counts = repmat(rlimVal-.1,[1 12]);
% % % % h = polarhistogram('BinEdges',edges,'BinCounts',counts,...
% % % %     'FaceColor','w','LineWidth',1,'FaceAlpha',1);
% % % % h.EdgeColor = 'none'; %accelColor;
% % % % hold on;

% shade for night
edges = linspace(-pi,pi,13);
counts = ones(1,12);
h = polarhistogram('BinEdges',edges,'BinCounts',counts,...
    'FaceColor','k','LineWidth',0.25,'FaceAlpha',0.7);
hold on;
h.EdgeColor = 'none';

% sun for day
edges = linspace(-pi,pi,size(sunYear,1)+1);
counts = sunYear(:,4) / 60 / 24;
h = polarhistogram('BinEdges',edges,'BinCounts',counts,...
    'FaceAlpha',1,'FaceColor',sunColor,'EdgeColor','none');

% temperature
counts = smooth(normalize(yearlyWeather,'range')+1);
edges = linspace(-pi,pi,numel(counts));
colors = parula(numel(counts));
colorLookup = linspace(1,2,size(colors,1));
for ii = 1:numel(counts)
    colorId = closest(colorLookup,counts(ii));
    polarplot([edges(ii) edges(ii)],[1 counts(ii)],'-','Color',colors(colorId,:),...
        'MarkerSize',15,'LineWidth',1);
    hold on;
end

iBehavior = 2; % accel
edges = linspace(-pi,pi,size(dailyMean,1)+1);
counts = normalize(inpaint_nans(dailyMean(:,iBehavior),4),'range')+1;
% % % % h = polarhistogram('BinEdges',edges,'BinCounts',counts,...
% % % %     'FaceColor',accelColor,'LineWidth',0.5,'FaceAlpha',0.5);
% % % % % h.DisplayStyle = 'stairs';
% % % % h.EdgeColor = 'none'; %accelColor;
% % % % hold on;
% outline accell
h = polarhistogram('BinEdges',edges,'BinCounts',counts,...
    'FaceColor',accelColor,'LineWidth',1,'FaceAlpha',0.5);
h.DisplayStyle = 'stairs';
% h.EdgeColor = accelColor;
h.EdgeColor = 'k';

h = polarplot(edges,ones(size(edges))*1.5); 
h.Color = repmat(0.75,[1,3]);
h.LineStyle = ':';

% black border at r=1
h = polarplot(edges,ones(size(edges))); 
h.Color = 'k';

pax = gca;
pax.ThetaZeroLocation = 'bottom';
pax.ThetaDir = 'clockwise';
pax.FontSize = 18;
pax.Layer = 'top';
% rlim([0 rlimVal]);
rticks([]);
pax.Color = [1 1 1];
% rticklabels({'','',''});
fontSize = 14;
% % % % text(-pi/2,1,'\leftarrow24 hrs','FontSize',fontSize,'Color','k','HorizontalAlignment','left');
text(0.67,1.5,'activity','FontSize',fontSize,'Color','k');
thetaticklabels(circshift(monthNames,6));
text(-pi/2-.01,1.5,['\leftarrow0',char(176),'C'],'FontSize',fontSize,'Color','k','HorizontalAlignment','left');
