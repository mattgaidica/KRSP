function addFigureLabels(h,varargin)
allAxes = findall(h,'type','axes');

labelMap = ['A';'B';'C';'D';'E';'F';'G'];
labelText = labelMap(1:numel(allAxes));

XY = repmat([-.25,1.1],[numel(allAxes),1]);
if nargin > 1
    XY = varargin{1};
end

pos = h.Position;
fontSize = 22;%round(pos(3)/50,2); % 16pts for 1200px figure

figure(h);
subplotIdx = flip(1:numel(allAxes)); % these are reversed
for iSubplot = 1:numel(allAxes)
    subplot(allAxes(subplotIdx(iSubplot)));
    text(XY(iSubplot,1),XY(iSubplot,2),labelText(iSubplot),'fontsize',fontSize,...
        'units','normalized','fontweight','bold');
end

