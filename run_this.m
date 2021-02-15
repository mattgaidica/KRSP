ff(1000,1000);
rows = 3;
cols = 3;

for iPlot = 1:rows*cols
    subplot(rows,cols,iPlot);
    plot(data(iPlot,:),'k');
end