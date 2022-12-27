function [nestFixed,shiftNest] = fixTempDelay(nest,odba,temp)
doPlot = 1;
if iscell(nest)
    nest = strcmp(nest,'Nest');
end

searchWindow = 60*10; % minutes

% % nestExits = find(diff(nest)==-1);
% % tempGrad = [];
% % kneePts = [];
% % jj = 0;
% % z_temp = normalize(temp);
% % for ii = 1:numel(nestExits)
% %     useRange = nestExits(ii)-searchWindow+1:nestExits(ii);
% %     if min(useRange) > 0
% %         if z_temp(useRange(1)) > 0 && minmax(z_temp(useRange)) > 1
% %             jj = jj + 1;
% %             tempGrad(jj,:) = temp(useRange);
% %             kneePts(jj) = find(gradient(temp(useRange)) < 0.2/60,1,'last');
% %         end
% %     end
% % end
% % meanGrad = mean(tempGrad);
% % t = linspace(-searchWindow,0,numel(meanGrad));
% % [~,kneeIdx] = knee_pt(t,meanGrad);
% % shiftNest = round(t(kneeIdx));
% % nestFixed = circshift(nest,shiftNest); % rotate
% % nestFixed(end+shiftNest+1:end) = repmat(nestFixed(end+shiftNest),[1,-shiftNest]); % fill in end
% % 
% % ff(300,400);
% % plot(t,tempGrad');
% % hold on;
% % plot(t,meanGrad,'k','linewidth',3);
% % xline(shiftNest,'r-','linewidth',2);
% % grid on;
% % xlabel('Time (s)');
% % ylabel('Temp (C)');
% % title(sprintf("Avg Delay = %i seconds",shiftNest));
% % fprintf("\nShift Nest: %i seconds\n\n",shiftNest);
% % 
% % return;

% ODBA method
searchIncs = [60,10,1];
searchRange = 1:searchIncs(1):searchWindow; % init
fprintf("Searching for minimum, backstepping %1.0f minutes\n",searchWindow/60);
for iSearch = 1:numel(searchIncs)
    fprintf("Narrowing: %i-%i seconds (lag)...\n",searchRange(1),searchRange(end));
    meanODBA = NaN(numel(searchRange),1);
    jj = 0;
    for ii = searchRange
        useRange = 1:numel(nest)-searchWindow;
        nestShift = circshift(nest,-ii);
        nestSelect = nestShift(useRange);
        odbaSelect = odba(useRange);
        jj = jj + 1;
        meanODBA(jj) = mean(odbaSelect(nestSelect==1));
    end
    [~,k] = min(meanODBA);
    if iSearch < numel(searchIncs)
        searchRange = searchRange(max([k-1,1])):searchIncs(iSearch+1):searchRange(min([k+1,jj]));
    end
end

shiftNest = -searchRange(k); % use in circshift
fprintf("\nShift Nest: %i seconds\n\n",shiftNest);
nestFixed = circshift(nest,shiftNest); % rotate
nestFixed(end+shiftNest+1:end) = repmat(nestFixed(end+shiftNest),[1,-shiftNest]); % fill in end
    
if doPlot
    ff(500,225);
    plot(searchRange,meanODBA,'k-','linewidth',2);
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

% this works, but performing gradient over chunky data is inaccurate
% it really needs to focus on the transition area!
% % searchRange = 1:searchIncs(1):searchWindow; % init
% % fprintf("Searching for minimum, backstepping %1.0f minutes\n",searchWindow/60);
% % for iSearch = 1:numel(searchIncs)
% %     fprintf("Narrowing: %i-%i seconds (lag)...\n",searchRange(1),searchRange(end));
% %     meanArr = NaN(numel(searchRange),1);
% %     jj = 0;
% %     for ii = searchRange
% %         useRange = 1:numel(nest)-searchWindow;
% %         nestShift = circshift(nest,-ii);
% %         nestSelect = nestShift(useRange);
% %         tempSelect = temp(useRange);
% %         gradRange = gradient(tempSelect(nestSelect==0)); % out of nest
% %         gradRange(sign(gradRange) >= 0) = NaN; % set non-negative gradient = NaN
% %         jj = jj + 1;
% %         meanArr(jj) = nanmean(gradRange);
% %     end
% %     [~,k] = min(meanArr);
% %     if iSearch < numel(searchIncs)
% %         searchRange = searchRange(max([k-1,1])):searchIncs(iSearch+1):searchRange(min([k+1,jj]));
% %     end
% % end
% % 
% % shiftNest = -searchRange(k); % use in circshift
% % fprintf("\nShift Nest: %i seconds\n\n",shiftNest);
% % nestFixed = circshift(nest,shiftNest); % rotate
% % nestFixed(end+shiftNest+1:end) = repmat(nestFixed(end+shiftNest),[1,-shiftNest]); % fill in end
% %     
% % if doPlot
% %     ff(500,225);
% %     plot(searchRange,meanArr,'k-','linewidth',2);
% %     hold on;
% %     xline(-shiftNest,'r','linewidth',2);
% %     set(gca,'fontsize',14);
% %     xlabel('Lag (sec)');
% %     grid on;
% %     title(sprintf('Minimizing Temp Gradient In Nest\nTemp lags Nest by %1.0f sec (%1.2f min)',...
% %         -shiftNest,-shiftNest/60));
% %     ylabel('ODBA');
% %     xlim([min(searchRange) max(searchRange)]);
% %     % saveas(gcf,fullfile(savePath,'nest-axy-lag.jpg'));
% % end