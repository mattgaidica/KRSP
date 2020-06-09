function compile_odba(filespath)
decimateBy = 100;

files = dir2(filespath,'*.csv','-r');
for iFile = 1:numel(files)
    readFile = fullfile(filespath,files(iFile).name);
    saveFile = [readFile,'.mat'];
    if ~exist(saveFile)
        disp(['working on ',readFile]);
        inputTable = readtable(readFile);
        nRange = 1:decimateBy:size(inputTable,1);
        T = table;
        if ismember('datetime',inputTable.Properties.VariableNames)
            dtData = inputTable.datetime(nRange);
        else
            dtData = inputTable.dtime(nRange);
        end
        if isa(dtData(1),'datetime')
            T.datetime = dtData;
        else
            T.datetime = datetime(dtData,'InputFormat','dd-MM-yyyy HH:mm:ss.00000');
        end
        T.odba = inputTable.odba(nRange);
        T.temp = inputTable.tempC(nRange);
        save(saveFile,'T');
    else
        disp([readFile,' already exists, skipping']);
    end
end