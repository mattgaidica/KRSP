xlFiles = dir2('/Volumes/GAIDICASSD/KRSP/KRSP Axy Data','-r','*.xlsx');
% remove hidden files
rmIds = [];
for ii = 1:size(xlFiles)
    [~,name,~] = fileparts(xlFiles(ii).name);
    if name(1) == '.'
        rmIds = [rmIds ii]; %#ok<AGROW> 
    end
end
xlFiles(rmIds,:) = [];
%%
foundAxyIds = [];
for ii = 1:size(xlFiles)
    csvFiles = dir2(xlFiles(ii).folder,'-r','*.csv');
    for jj = 1:numel(csvFiles)
        [~,name,ext] = fileparts(csvFiles(jj).name);
        foundAxyIds = [foundAxyIds;find(strcmp(T_AxyFiles.filename,[name,ext]))]; %#ok<AGROW> 
    end
end
foundAxyIds = unique(foundAxyIds);
clc
for ii = 1:size(T_AxyFiles,1)
    if ismember(ii,foundAxyIds)
        continue;
    end
    fprintf("%s\n",T_AxyFiles.folder(ii));
end