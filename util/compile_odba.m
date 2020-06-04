function compile_odba(filespath)
decimateBy = 100;

files = dir2(filespath,'*.csv','-r');
for iFile = 1:numel(files)
    readFile = fullfile(files(iFile).folder,files(iFile).name);
    saveFile = [readFile,'.mat'];
    if ~exist(saveFile)
        disp(['working on ',readFile]);
        inputTable = readtable(readFile);
        nRange = 1:decimateBy:size(inputTable,1);
        T = table;
        if isa(inputTable.datetime(1),'datetime')
            T.datetime = inputTable.datetime(nRange);
        else
            T.datetime = datetime(inputTable.datetime(nRange),'InputFormat','dd-MM-yyyy HH:mm:ss.00000');
        end
        T.odba = inputTable.odba(nRange);
        T.temp = inputTable.tempC(nRange);
        save(saveFile,'T');
    else
        disp([readFile,' already exists, skipping']);
    end
end