axyDbFile = '/Volumes/GAIDICASSD/KRSP/KRSP Axy Data/AxyDatabase.xlsx';
axyPath = '/Users/matt/Documents/MATLAB/KRSP/Axy';
opts = detectImportOptions(axyDbFile);

% opts = setvartype(opts,opts.VariableNames(contains(opts.VariableNames,{'shake','time'})),'datetime');
T_AxyDB = readtable(axyDbFile,opts);
fprintf("total filenames: %i/%i rows\n",sum(strcmp(T_AxyDB.filename,"")==0),size(T_AxyDB,1));
save(fullfile(axyPath,'T_AxyDB.mat'),'T_AxyDB','axyDbFile','opts');

%% re-export axy files?
storagePath = '/Volumes/GAIDICASSD/KRSP';
clc
fullfiles = string;
for ii = 1:size(T_AxyDB,1)
    if strcmp(T_AxyDB.filename{ii},"") == 0
        csvFiles = dir3(fullfile(storagePath,T_AxyDB.folder{ii}),T_AxyDB.filename{ii});
        if size(csvFiles,1) == 1
            fullfiles(ii,:) = string(fullfile(csvFiles.folder,T_AxyDB.filename{ii}));
        else
            fullfiles(ii,:) = "";
            fprintf("error %i %s\n",ii,T_AxyDB.filename{ii});
        end
    else
        fullfiles(ii,:) = "";
    end
end


%% find squirrel_id
% krsp_squirrel_file = '/Volumes/GAIDICASSD/KRSP/krsp_squirrel.csv';
% opts = detectImportOptions(krsp_squirrel_file);
% opts = setvartype(opts,opts.VariableNames(contains(opts.VariableNames,{'tag'})),'string');
% krsp_squirrel = readtable(krsp_squirrel_file,opts);
