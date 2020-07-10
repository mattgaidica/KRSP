function Tawake = make_Tawake(T)
Tawake = table;
nId = T.awake == 1;
diff_nId = diff(nId);
changes = [1;find(diff_nId)+1];
Tawake.datetime = T.datetime(changes);
Tawake.awake = T.awake(changes);
% for iChange = 1:numel(changes)
%     startId = changes(iChange);
%     if iChange == numel(changes)
%         endId = size(inputTable,1);
%     else
%         endId = changes(iChange+1);
%     end
%     Tstat.odba_max(iChange) = max(inputTable.odba(startId:endId));
%     Tstat.odba_mean(iChange) = mean(inputTable.odba(startId:endId));
%     Tstat.odba_sum(iChange) = sum(inputTable.odba(startId:endId));
%     Tstat.odba_med(iChange) = median(inputTable.odba(startId:endId));
%     Tstat.odba_std(iChange) = std(inputTable.odba(startId:endId));
%     Tstat.temp_mean(iChange) = mean(inputTable.tempC(startId:endId));
% end