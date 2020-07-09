function compile_odba(filespath,qualifier)

warning('off','all');
files = dir2(filespath,['*',qualifier,'.csv'],'-r');
for iFile = 1:numel(files)
    readFile = fullfile(filespath,files(iFile).name);
    disp(['working on: ...',readFile(end-40:end)]);
    inputTable = readtable(readFile);
    T = table;
    if ismember('datetime',inputTable.Properties.VariableNames)
        dtData = inputTable.datetime;
    else
        dtData = inputTable.dtime;
    end
    if ~isa(dtData(1),'datetime')
        dtData = datetime(dtData,'InputFormat','dd-MM-yyyy HH:mm:ss.00000');
    end
    medDays = median(diff(unique(day(dtData,'dayofyear'))));
    if ~isnan(medDays) && medDays ~= 1 % NaN is single day && days progress by 1
        dtData = datetime(year(dtData),day(dtData),month(dtData),...
            hour(dtData),minute(dtData),second(dtData));
        fprintf('formatting... ');
    end
    fprintf('%3.0f days recorded\n',days(dtData(end)-dtData(1)));
    
    Fs = 1 / seconds(median(diff(dtData))); % 1 / period
    decimateBy = 60*Fs; % compress to 1 minute
    nRange = 1:decimateBy:size(inputTable,1);
    T.datetime = dtData(nRange);
    T.odba = inputTable.odba(nRange);
    odba_max = movmax(inputTable.odba,decimateBy);
    T.odba_max = odba_max(nRange);
    
    if ismember('tempC',inputTable.Properties.VariableNames)
        T.temp = inputTable.tempC(nRange);
        T.nest = inputTable.Nest2(nRange);
    end
    
    Tstat = table;
    nId = strcmp(inputTable.Nest2,'Nest');
    diff_nId = diff(nId);
    changes = [1;find(diff_nId)+1];
    Tstat.datetime = dtData(changes);
    Tstat.nest = inputTable.Nest2(changes);
    for iChange = 1:numel(changes)
        startId = changes(iChange);
        if iChange == numel(changes)
            endId = size(inputTable,1);
        else
            endId = changes(iChange+1);
        end
        Tstat.odba_max(iChange) = max(inputTable.odba(startId:endId));
        Tstat.odba_mean(iChange) = mean(inputTable.odba(startId:endId));
        Tstat.odba_sum(iChange) = sum(inputTable.odba(startId:endId));
        Tstat.odba_med(iChange) = median(inputTable.odba(startId:endId));
        Tstat.odba_std(iChange) = std(inputTable.odba(startId:endId));
        Tstat.temp_mean(iChange) = mean(inputTable.tempC(startId:endId));
    end
    
    saveFile = strrep(readFile,'.csv',['__',datestr(dtData(1),'yyyymmdd'),'.mat']);
    disp(['saving ',saveFile(end-40:end)]);
    save(saveFile,'T','Tstat');
end
warning('on','all');