function [shiftNest,nestFixed] = findTempDelay(nest,odba)
doPlot = 1;
if iscell(nest)
    nest = strcmp(nest,'Nest');
end

tic;
searchWindow = 60*10; % minutes
searchIncs = [60,10,1];
searchRange = 1:searchIncs(1):searchWindow; % init
fprintf("Searching for minimum, backstepping %1.0f minutes\n",searchWindow/60);
for iSearch = 1:numel(searchIncs)
    fprintf("Narrowing: %i-%i seconds (lag)...\n",searchRange(1),searchRange(end));
    sumODBA = NaN(numel(searchRange),1);
    jj = 0;
    for ii = searchRange
        useRange = 1:numel(nest)-searchWindow;
        nestShift = circshift(nest,-ii);
        nestSelect = nestShift(useRange);
        odbaSelect = odba(useRange);
        jj = jj + 1;
        sumODBA(jj) = mean(odbaSelect(nestSelect==1));
    end
    disp(jj);
    [~,k] = min(sumODBA);
    if iSearch < numel(searchIncs)
        searchRange = searchRange(max([k-1,1])):searchIncs(iSearch+1):searchRange(min([k+1,jj]));
    end
end
toc

% % % % sumODBA = NaN(searchWindow,1);
% % % % for ii = 1:searchWindow
% % % %     if ii == 1
% % % %         fprintf("Searching for minimum, backstepping %1.0f minutes\n0%%",searchWindow/60);
% % % %     elseif mod(ii,50) == 0
% % % %         fprintf("\n%i%%",round(100*ii/searchWindow));
% % % %     else
% % % %         fprintf(".");
% % % %     end
% % % %     useRange = 1:numel(nest)-searchWindow;
% % % %     nestShift = circshift(nest,-ii);
% % % %     nestSelect = nestShift(useRange);
% % % %     odbaSelect = odba(useRange);
% % % %     sumODBA(ii) = mean(odbaSelect(nestSelect==1));
% % % % end
% % % % [v,k] = min(sumODBA);

shiftNest = -searchRange(k); % use in circshift
fprintf("\nShift Nest: %i seconds\n\n",shiftNest);
nestFixed = circshift(nest,shiftNest); % rotate
nestFixed(end+shiftNest+1:end) = repmat(nestFixed(end+shiftNest),[1,-shiftNest]); % fill in end
    
if doPlot
    close all
    ff(500,225);
    plot(searchRange,sumODBA,'k-','linewidth',2);
    hold on;
    xline(-shiftNest,'r','linewidth',2);
    set(gca,'fontsize',14);
    xlabel('Lag (sec)');
    grid on;
    title(sprintf('Minimizing ODBA In Nest\nODBA lags Nest by %1.0f sec (%1.2f min)',...
        -shiftNest,-shiftNest/60));
    ylabel('ODBA');
    xlim([min(searchRange) max(searchRange)]);
    % saveas(gcf,fullfile(savePath,'nest-axy-lag.jpg'));
end