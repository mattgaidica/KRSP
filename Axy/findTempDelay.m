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
for ii = 1:numel(enterLocs)
    fprintf("Location %i/%i\n",ii,numel(enterLocs));
    exitIdx = find(exitLocs > enterLocs(ii),1,'first');
    if isempty(exitIdx)
        break; % done
    end
    useRange = enterLocs(ii):exitLocs(exitIdx);
    if numel(useRange) < 60
        continue; % too short
    end
    meanArr = [];
    for jj = 0:searchWindow-1
        meanArr(jj+1) = mean(odba(useRange-jj));
    end
    [~,k] = min(meanArr);
    minCount = minCount + 1;
    minArr(minCount) = k;
end
shiftNest = median(minArr);
nestFixed = circshift(strcmp(T.Nest,'Nest'),-shiftNest); % rotate
nestFixed(end-shiftNest+1:end) = repmat(nestFixed(end-shiftNest),[1,shiftNest]); % fill in end
fprintf("\nShift Nest: %i\n\n",shiftNest);

% old method
% % % % sumODBA = [];
% % % % for ii = 1:searchWindow
% % % %     if ii == 1
% % % %         fprintf("Searching for minimum, backstepping %1.0f minutes\n0%%",searchWindow/60);
% % % %     elseif mod(ii,50) == 0
% % % %         fprintf("\n%i%%",round(100*ii/searchWindow));
% % % %     else
% % % %         fprintf(".");
% % % %     end
% % % %     useRange = 1:numel(nestFixed)-searchWindow;
% % % %     nestShift = circshift(nestFixed,-ii);
% % % %     nestSelect = nestShift(useRange);
% % % %     odbaSelect = odba(useRange);
% % % %     sumODBA(ii) = mean(odbaSelect(nestSelect==1));
% % % % end
% % % % [v,k] = min(sumODBA);
% % % % shiftNest = -k; % use in circshift
% % % % nestFixed = circshift(strcmp(T.Nest,'Nest'),shiftNest); % rotate
% % % % if k > 0 % positive k
% % % %     nestFixed(end-k+1:end) = repmat(nestFixed(end-k),[1,k]); % fill in end
% % % % else % negative k
% % % %     nestFixed(1:-k) = repmat(nestFixed(-k),[1,-k]);
% % % % end

if doPlot
    t = linspace(0,searchWindow,searchWindow);
    close all
    ff(500,225);
    plot(t,sumODBA,'k-','linewidth',2);
    hold on;
    plot(t(k),v,'x','markerSize',20);
    set(gca,'fontsize',14);
    xlabel('Lag (sec)');
    grid on;
    title(sprintf('Minimizing ODBA In Nest\nODBA lags Nest by %1.0f sec (%1.2f min)',k,k/60));
    ylabel('ODBA');
    % saveas(gcf,fullfile(savePath,'nest-axy-lag.jpg'));
end