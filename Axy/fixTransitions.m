function binNestSense = fixTransitions(binNestSense,sense)
doPlot = false;
maxShift = 60*5; % window
transLocs = [1;find(abs(diff(binNestSense))==1);numel(binNestSense)]; % pad to rm logic below
% the norm pdf weighs location nearest to center as more important,
% could also consider using relative heights?
weightDist = normalize(normpdf(linspace(-2,2,maxShift*2)),'range',[0 1]);
for ii = 2:numel(transLocs)-1 % exclude pads
    % the following is to make sure numel(useRange) is maintained but
    % data from non-existant ids OR nearby trasitions is ignored
    useRange = transLocs(ii)-maxShift+1:transLocs(ii)+maxShift;
    usefulStart = max([1,transLocs(ii-1)]);
    usefulEnd = min([numel(binNestSense),transLocs(ii+1)]);
    fillStart = find(useRange < usefulStart,1,'last');
    fillEnd = find(useRange >= usefulEnd,1,'first');
    if isempty(fillStart); fillStart = 1; end
    if isempty(fillEnd); fillEnd = numel(useRange); end
    if fillEnd-fillStart < 120; continue; end
    % data used for sensing
    nanData = NaN(size(useRange));
    use_tempZ = nanData;
    use_tempZ(fillStart:fillEnd) = sense.tempZ(useRange(fillStart):useRange(fillEnd));
    % use gradient instead of diff to maintain dimensions
    sensorData1 = normalize(weightDist .* smoothdata(gradient(use_tempZ),'gaussian',60));
    if binNestSense(transLocs(ii)) == 1 % in to out, search for neg peak
        sensorData1 = -sensorData1;
    end
    sensorData2 = normalize(gradient(sensorData1));
    [~,k] = max(sensorData2); % find peak
    % fill data
    newData = binNestSense(useRange);
    if k >= maxShift
        newData(maxShift:k) = binNestSense(transLocs(ii));
    else
        newData(k:maxShift) = ~binNestSense(transLocs(ii));
    end
    
    if doPlot
        close all
        ff(600,500);
        yyaxis left;
        plot(sensorData2,'k-','linewidth',2);
        hold on;
        plot(k,sensorData2(k),'k.','markerSize',20);
        plot(sensorData1,'k:');
        
        yyaxis right;
        plot(use_tempZ,'r-','linewidth',2);
        hold on;
        plot(binNestSense(useRange),'g-','linewidth',1);
        plot(newData+1.1,'g-','linewidth',3);
        plot(weightDist,'b-');
        xlim([1 numel(useRange)]);
        grid on;
        hold off;
    end
    binNestSense(useRange) = newData; % overwrite
end