function compile_odba(filespath)
decimateBy = 100;

files = dir2(filespath,'*.csv','-r');
for iFile = 1:numel(files)
    readFile = fullfile(filespath,files(iFile).name);
    disp(['working on ',readFile]);
    inputTable = readtable(readFile);
    nRange = 1:decimateBy:size(inputTable,1);
    T = table;
    if ismember('datetime',inputTable.Properties.VariableNames)
        dtData = inputTable.datetime(nRange);
    else
        dtData = inputTable.dtime(nRange);
    end
    if days(dtData(end)-dtData(1)) > 90 % assume longest deploy
        disp('...date converted');
        dtData = datetime(year(dtData),day(dtData),month(dtData),...
            hour(dtData),minute(dtData),second(dtData));
    end
    if isa(dtData(1),'datetime')
        T.datetime = dtData;
    else
        T.datetime = datetime(dtData,'InputFormat','dd-mmm-yyyy HH:MM:SS');
    end
    T.odba = inputTable.odba(nRange);
    T.temp = inputTable.tempC(nRange);
    % ADD FIRST DATETIME to filename
    saveFile = [readFile,'_',datestr(dtData(1),'yyyymmdd'),'.mat'];
    save(saveFile,'T');
end