function cb = cbAside(ax,labelText,cbColor,caxisVals)
cbOffset = .01;
cb = colorbar(ax,'Location','east');
% cb.FontSize = 8;
gcaPos = ax.Position;
cbPos = cb.Position;
set(cb,'position',[gcaPos(1) + gcaPos(3) + cbOffset gcaPos(2) cbPos(3)/2 gcaPos(4)]);
set(cb,'YAxisLocation','right');
cb.Color = 'r';

if exist('labelText')
    cb.Label.String = labelText;
end

if exist('cbColor')
    cb.Color = cbColor;
else
    cb.Color = 'r';
end

if exist('caxisVals')
    caxis(ax,caxisVals);
else
    if sum(abs(caxis)) > 10
        caxis(ax,round(caxis));
    end
end
cb.Ticks = caxis;
cb.TickLabels = compose("%1.1e",caxis);
cb.Label.Position(1) = cbPos(1)+0.1;

cb.Box = 'off';