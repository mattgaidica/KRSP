axyDbFile = '/Volumes/GAIDICASSD/KRSP/KRSP Axy Data/AxyDatabase.xlsx';
opts = detectImportOptions(axyDbFile);

opts = setvartype(opts,opts.VariableNames(contains(opts.VariableNames,{'shake','time'})),'datetime');
T_AxyDB = readtable(axyDbFile,opts);
fprintf("total filenames: %i/%i rows\n",sum(strcmp(T_AxyDB.filename,"")==0),size(T_AxyDB,1));
save('T_AxyDB.mat','T_AxyDB','storagePath','axyDbFile','opts');

%% re-export axy files
storagePath = '/Volumes/GAIDICASSD/KRSP';
clc
fullfiles = string;
for ii = 1:size(T_AxyDB,1)
    if strcmp(T_AxyDB.filename{ii},"") == 0
        csvFiles = dir2(fullfile(storagePath,T_AxyDB.root_dir{ii}),'-r',T_AxyDB.filename{ii});
        if size(csvFiles,1) == 1
            fullfiles(ii,:) = string(fullfile(csvFiles.folder,T_AxyDB.filename{ii}));
        else
            fullfiles(ii,:) = "";
            fprintf("error %i %s\n",ii,T_AxyDB.filename{ii});
        end
    end
end


%% find squirrel_id
% krsp_squirrel_file = '/Volumes/GAIDICASSD/KRSP/krsp_squirrel.csv';
% opts = detectImportOptions(krsp_squirrel_file);
% opts = setvartype(opts,opts.VariableNames(contains(opts.VariableNames,{'tag'})),'string');
% krsp_squirrel = readtable(krsp_squirrel_file,opts);
clc
filenames = string;
for ii = 1:size(T_AxyDB,1)
    if strcmp(T_AxyDB.filename{ii},"") == 0
        matchIdx = find(strcmp(T_AxyDB.filename{ii},T_AxyFiles.filename));
        if isempty(matchIdx)
            fprintf("missing: %s (%i)\n",T_AxyDB.filename{ii},ii+1);
        end
        if numel(matchIdx) > 1
            fprintf("multiple: %s (%i)\n",T_AxyDB.filename{ii},ii+1);
        end
    end
end
%%
csvFiles = dir2('/Volumes/GAIDICASSD/KRSP/KRSP Axy Data/Dantzer Lab Axy Data/AXY Data 2017','-r','*.csv');
filenames = "";
middens = "";
axys = "";
for ii = 1:size(csvFiles)
    [path,name,ext] = fileparts(csvFiles(ii).name);
    nameSplit = strsplit(name,"_");
    axys(ii,:) = string(nameSplit{1}(4:end));
    splitPath = strsplit(path);
    middens(ii,:) = string(splitPath{1});
    filenames(ii,:) = string([name,ext]);
end
%%
storagePath = '/Volumes/GAIDICASSD/KRSP';
csvFiles = dir2(storagePath,'-r','*.csv');
rmIds = [];
for ii = 1:size(csvFiles)
    [~,name] = fileparts(csvFiles(ii).name);
    if strcmp(name(1),'.')
        rmIds = [rmIds ii];
    end
end
csvFiles(rmIds) = [];

%%
clc
for ii = 1:size(T_AxyDB,1)
    matchIdx = find(contains({csvFiles.name},T_AxyDB.filename{ii}) & contains({csvFiles.name},T_AxyDB.root_dir{ii}));
    if numel(matchIdx) > 1 && strcmp(T_AxyDB.filename{ii},"") == 0
        for jj = 1:numel(matchIdx)
            fprintf("%s - %s\n",T_AxyDB.filename{ii},T_AxyDB.root_dir{ii});
        end
        fprintf("\n");
    end
end