function zOffset = findZOffset(temp,odba,nest)
doPlot = true;

% this fixes the nest guess from temp only then minimizes ODBA in nest
% while maximizing ODBA out of nest
rmShortNest = removeShortTransitions(nest,60); % optional
shiftNest = fixTempDelay(rmShortNest,odba,temp); % re-align nest
% most odba in-nest should be Z<0, most out Z>0
addVals = linspace(-1,1,100);
resArr = NaN(size(addVals));

for ii = 1:numel(addVals)
    newComp = normalize(odba) + addVals(ii);
    resArr(ii,1) = sum(sign(newComp(shiftNest==0))==1)/sum(shiftNest==0);
    resArr(ii,2) = sum(sign(newComp(shiftNest==1))==-1)/sum(shiftNest==1);
end
nSmooth = round(numel(addVals)*.1);
sm1 = smoothdata(resArr(1:ii,1),'gaussian',nSmooth);
sm2 = smoothdata(resArr(1:ii,2),'gaussian',nSmooth);
addIdx = find(sm1>sm2,1,'first');
if isempty(addIdx)
    zOffset = 0;
else
    zOffset = addVals(addIdx);
end


if doPlot
    colors = lines(5);
    ff(500,300);
    plot(addVals(1:ii),resArr(1:ii,1),'k:','linewidth',1);
    hold on;
    ln1 = plot(addVals(1:ii),sm1,'k-','linewidth',2);
    plot(addVals(1:ii),resArr(1:ii,2),':','color',colors(5,:),'linewidth',1);
    ln2 = plot(addVals(1:ii),sm2,'-','color',colors(5,:),'linewidth',2);
    grid on;
    if ~isempty(addIdx)
        xline(addVals(addIdx),'r-','linewidth',2);
    end

    title(sprintf('max out-nest/min in-nest Z=%1.2f',zOffset));
    xlabel('Z-offset');
    ylabel('ODBA');
    legend([ln1,ln2],["In-nest ODBA","Out-nest ODBA"]);
    set(gca,'fontsize',14);
end