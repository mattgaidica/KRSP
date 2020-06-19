function compile_odba(filespath)
decimateBy = 60;

warning('off','all');
files = dir2(filespath,'*.csv','-r');
for iFile = 1:numel(files)
    readFile = fullfile(filespath,files(iFile).name);
    disp(['working on: ...',readFile(end-40:end)]);
    inputTable = readtable(readFile);
    nRange = 1:decimateBy:size(inputTable,1);
    T = table;
    if ismember('datetime',inputTable.Properties.VariableNames)
        dtData = inputTable.datetime(nRange);
    else
        dtData = inputTable.dtime(nRange);
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
    T.datetime = dtData;
    T.odba = inputTable.odba(nRange);
    T.temp = inputTable.tempC(nRange);
    T.Nest = inputTable.Nest2(nRange);
    saveFile = strrep(readFile,'.csv',['__',datestr(dtData(1),'yyyymmdd'),'.mat']);
    disp(['saving ',saveFile(end-40:end)]);
    save(saveFile,'T');
    warning('on','all');
end