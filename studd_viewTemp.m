% dataFile = '/Users/matt/Documents/Data/KRSP/doi_10.5061_dryad.1s1m8r7__v1/squirrelAxy_decisionTree.csv';
% A = readtable(dataFile);
monthNames = {'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
[v,k] = sort(A.TIME);
% [C,ia,ic] = studd_help_labelSq(A);

close all     

for ii = 1:numel(C)
    theseTemps = A.tempC(ic==ii);
    theseTimes = A.TIME(ic==ii);
    theseBeh = A.BEHAV(ic==ii);
%     [~,k] = sort(A.TIME(ic==ii));
%     theseData = theseData(k); % sort by time
    figure;
    daySpan = min(day(theseTimes,'dayofyear')):max(day(theseTimes,'dayofyear'));
    colors = parula(numel(daySpan));
    dayCount = 0;
    for iDay = daySpan
        dayCount = dayCount + 1;
        todayTime = theseTimes(day(theseTimes,'dayofyear') == iDay);
        if ~isempty(todayTime)
            today_seconds = hour(todayTime)*60*60 + minute(todayTime)*60 + second(todayTime);
            todayTemp = theseTemps(day(theseTimes,'dayofyear') == iDay);
            todayBeh = theseBeh(day(theseTimes,'dayofyear') == iDay);
            [sortedSeconds,k] = sort(today_seconds);
            sortedTemps = todayTemp(k);
            sortedBeh = todayBeh(k);
            for jj = 1:numel(sortedTemps)
                marker = '.';
                if strcmp(sortedBeh{jj},'Nest')
                    marker = '*';
                end
                plot(sortedSeconds(jj)/(60*60),sortedTemps(jj),marker,'markerSize',10,'color',colors(dayCount,:));
                hold on;
            end
        end
    end
end

% theta = [];
% rho = [];
% used = [];
% dailyTemp = [];
% for ii = 1:numel(v) % each day
%     % find mean for that day
%     d = v(ii);
%     
%     if isempty(used) || ~ismember(day(d,'dayofyear'),day(used,'dayofyear'))% || ~ismember(year(d),year(used))
%         used = [used;d];
%         dailyTemp(numel(used)) = mean(A.tempC(ismember(day(v,'dayofyear'),day(d,'dayofyear'))));
%         theta(numel(used)) = (2*pi/300) * day(d,'dayofyear');
%     end
% end
% 
% figure;
% polarplot(theta,dailyTemp);
% 
% thetaticklabels(circshift(monthNames,6));
% pax = gca;
% pax.ThetaZeroLocation = 'bottom';
% pax.ThetaDir = 'clockwise';
% pax.FontSize = 18;
% pax.Layer = 'top';