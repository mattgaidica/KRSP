function axy = april_tempFilter(axy)
% remove large deviations due to low power mode
axy.temp = movmedian(axy.temp,10);

% % % % function axy = april_tempFilter(axy, posThreshold, negThreshold)
% % % % % Ensure the datetime is sorted
% % % % axy = sortrows(axy, 'datetime');
% % % % 
% % % % % Extract times and temperatures
% % % % datetime = axy.datetime;
% % % % temp = axy.temp;
% % % % 
% % % % % Find unique 10 second intervals
% % % % m = min(datetime):seconds(10):max(datetime);
% % % % m = m';
% % % % 
% % % % % Find closest index of original data to the 10 sec interval
% % % % [~, idx] = min(abs(datetime - m'), [], 2);
% % % % 
% % % % % Group temperature data according to these 10 second intervals
% % % % id = idx;
% % % % sumTable = table(id, temp);
% % % % 
% % % % % Calculate mean temperature and difference in temperature
% % % % sumTable = varfun(@mean, sumTable, 'GroupingVariables', 'id', 'InputVariables', 'temp');
% % % % sumTable.difftemp1 = [NaN; diff(sumTable.mean_temp)];
% % % % sumTable.difftemp2 = [diff(sumTable.mean_temp); NaN];
% % % % 
% % % % % Filter and correct temperature
% % % % mask = ((sumTable.difftemp1 >= posThreshold | sumTable.difftemp1 <= negThreshold) & ...
% % % %         (sumTable.difftemp2 >= posThreshold | sumTable.difftemp2 <= negThreshold));
% % % % corrected_temp = (circshift(sumTable.mean_temp, -1) + circshift(sumTable.mean_temp, 1)) / 2;
% % % % sumTable.tempC = sumTable.mean_temp;
% % % % sumTable.tempC(mask) = corrected_temp(mask);
% % % % 
% % % % % Merge the corrected temperature back with the original data
% % % % axy = join(axy, sumTable(:, {'id', 'tempC'}), 'Keys', 'id');
% % % % 
% % % % % Drop the intermediate id column
% % % % axy.id = [];
