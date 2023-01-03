function [axyLogIds,no_axyLogIds,logFileList] = searchForLogs(T_AxyFiles,rootDir)
xlFiles = dir2(rootDir,'-r','*.xlsx');
% remove hidden files
rmIds = [];
for ii = 1:size(xlFiles)
    [~,name,~] = fileparts(xlFiles(ii).name);
    if name(1) == '.'
        rmIds = [rmIds ii]; %#ok<AGROW>
    end
end
xlFiles(rmIds,:) = [];

axyLogIds = [];
logFileList = string;
fCount = 0;
for ii = 1:size(xlFiles)
    csvFiles = dir2(xlFiles(ii).folder,'-r','*.csv');
    [~,name,ext] = fileparts(xlFiles(ii).name);
    logFilename = [name,ext];
    for jj = 1:numel(csvFiles)
        [~,name,ext] = fileparts(csvFiles(jj).name);
        axyId = find(strcmp(T_AxyFiles.filename,[name,ext]));
        for kk = 1:numel(axyId)
            fCount = fCount + 1;
            axyLogIds(fCount) = axyId(kk); %#ok<AGROW>
            logFileList(fCount,:) = string(fullfile(xlFiles(ii).folder,logFilename));
        end
    end
end
clc
no_axyLogIds = [];
for ii = 1:size(T_AxyFiles,1)
    if ismember(ii,axyLogIds)
        continue;
    end
    no_axyLogIds = [no_axyLogIds;ii];
    fprintf("No log: %s\n",fullfile(T_AxyFiles.folder(ii),T_AxyFiles.filename(ii)));
end