function categorize_odba(filespath,qualifier)
strFilts = {'Nest','Out'};
files = dir(fullfile(filespath,['*',qualifier,'.csv'])); % only this dir
% files = dir2(filespath,['*',qualifier,'.csv'],'-r');
warning('off','all');
nSmooth = 91; % !! DEPENDS ON FS - from Studd 2019
for iFile = 1:numel(files)
    readFile = fullfile(filespath,files(iFile).name);
    disp(readFile(end-60:end));
    inputTable = readtable(readFile);
    if isempty(inputTable) || ~ismember('Var1',inputTable.Properties.VariableNames)
        disp('table empty!');
        continue;
    end
    % raw data
    if ismember('Var1',inputTable.Properties.VariableNames) &&...
            ~ismember('odba',inputTable.Properties.VariableNames)
        inputTable.datetime = inputTable.Var1;
        X = inputTable.Var2;
        Y = inputTable.Var3;
        Z = inputTable.Var4;
        inputTable.odba = abs(X-medfilt1(X,nSmooth))...
            + abs(Y-medfilt1(Y,nSmooth))...
            + abs(Z-medfilt1(Z,nSmooth));
        % when is this condition used?
        tempDiff = max(inputTable.Var5) - min(inputTable.Var5);
        if std(inputTable.Var5) > 2 && tempDiff > 5 && tempDiff < 90
            inputTable.temp = inputTable.Var5;
        elseif ismember('Var6',inputTable.Properties.VariableNames)
            inputTable.temp = inputTable.Var6;
        else
            error('can not find temp data');
        end
    end
    % find each day, apply kmeans
    [doys,IA,IC] = unique(day(inputTable.datetime,'dayofyear'));
    disp(['no. days: ',num2str(numel(doys))]);
    nearestDayIds = [1;IA(3:end-1);numel(IC)]; % never less than a full day of data
    outputTable = table;
    outputTable.datetime = inputTable.datetime;
    outputTable.tempC = inputTable.temp;
    outputTable.odba = inputTable.odba;
    for nId = 1:numel(nearestDayIds)-1
        nRange = nearestDayIds(nId):nearestDayIds(nId+1);
        thisData = inputTable.temp(nRange);
        idx = kmeans(thisData,2); % 1==in, 2==out !! HOW do I know?
        if mean(thisData(idx==1)) > mean(thisData(idx==2))
            outputTable.Nest2(nRange(idx==1)) = strFilts(1);
            outputTable.Nest2(nRange(idx==2)) = strFilts(2);
        else
            outputTable.Nest2(nRange(idx==2)) = strFilts(1);
            outputTable.Nest2(nRange(idx==1)) = strFilts(2);
        end
    end
    if ~isfolder(fullfile(filespath,'nest'))
        mkdir(fullfile(filespath,'nest'))
    end
    writeFile = fullfile(filespath,'nest',strrep(files(iFile).name,'.csv','_nest.csv'));
    writetable(outputTable,writeFile);
end
warning('on','all');