function zOffset = findZOffset(temp,odba,nest)
% this fixes the nest guess from temp only then minimizes ODBA in nest
% while maximizing ODBA out of nest
rmShortNest = removeShortTransitions(nest,60); % optional
shiftNest = fixTempDelay(rmShortNest,odba,temp); % re-align nest
% most odba in-nest should be Z<0, most out Z>0
addVals = linspace(0,1,100);
resArr = NaN(size(addVals));
addIdx = 1;
for ii = 1:numel(addVals)
    newComp = normalize(odba) + addVals(ii);
    resArr(ii,1) = sum(sign(newComp(shiftNest==0))==1)/sum(shiftNest==0);
    resArr(ii,2) = sum(sign(newComp(shiftNest==1))==-1)/sum(shiftNest==1);
    if resArr(ii,1) >= resArr(ii,2)
        addIdx = ii;
        break;
    end
end
ff(300,200);
plot(addVals(1:ii),resArr(1:ii,1));
hold on;
plot(addVals(1:ii),resArr(1:ii,2));
grid on;

zOffset = addVals(addIdx);
title(sprintf('max out-nest/min in-nest Z=%1.2f',zOffset));
xlabel('Z-offset');
ylabel('ODBA');