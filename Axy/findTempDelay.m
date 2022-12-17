% function [shiftNest,nestFixed] = findTempDelay(T)
doPlot = 1;

searchWindow = 60*10; % minutes
nestFixed = strcmp(T.Nest,'Nest');
odba = T.odba;
diff_nest = diff(nestFixed);
enterLocs = find(diff_nest==1);
exitLocs = find(diff_nest==-1);

minArr = [];
minCount = 0;
ff(600,400);
for ii = 1:numel(enterLocs)
    fprintf("Location %i/%i\n",ii,numel(enterLocs));
    exitIdx = find(exitLocs > enterLocs(ii),1,'first');
    if isempty(exitIdx)
        break; % done
    end
    useRange = enterLocs(ii):exitLocs(exitIdx);
    if numel(useRange) < searchWindow
        continue; % too short
    end
    meanArr = [];
    for jj = 0:searchWindow-1
        meanArr(jj+1) = mean(odba(useRange-jj));
    end
    plot(meanArr);
    hold on;
    [~,k] = min(meanArr);
    minCount = minCount + 1;
    minArr(minCount) = k;
end
shiftNest = median(minArr);
nestFixed = circshift(strcmp(T.Nest,'Nest'),-shiftNest); % rotate
nestFixed(end-shiftNest+1:end) = repmat(nestFixed(end-shiftNest),[1,shiftNest]); % fill in end
fprintf("\nShift Nest: %i\n\n",shiftNest);

if doPlot
    ff(500,400);
    plot(sort(minArr),'k-','linewidth',2);
    
end