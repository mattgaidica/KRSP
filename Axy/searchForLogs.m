function [axyLogIds,no_axyLogIds] = searchForLogs(T_AxyFiles,rootDir)
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
for ii = 1:size(xlFiles)
    csvFiles = dir2(xlFiles(ii).folder,'-r','*.csv');
    for jj = 1:numel(csvFiles)
        [~,name,ext] = fileparts(csvFiles(jj).name);
        axyLogIds = [axyLogIds;find(strcmp(T_AxyFiles.filename,[name,ext]))]; %#ok<AGROW>
    end
end
axyLogIds = unique(axyLogIds);
clc
no_axyLogIds = [];
for ii = 1:size(T_AxyFiles,1)
    if ismember(ii,axyLogIds)
        continue;
    end
    no_axyLogIds = [no_axyLogIds;ii];
    fprintf("No log: %s\n",fullfile(T_AxyFiles.folder(ii),T_AxyFiles.filename(ii)));
end