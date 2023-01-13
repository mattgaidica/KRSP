% load
storagePath = '/Volumes/GAIDICASSD/KRSP/Axy Database';
axyDbFile = 'AxyDatabase.xlsx';
opts = detectImportOptions(fullfile(storagePath,axyDbFile));

% override opts
variableNames = opts.VariableNames;
[variableTypes{1:numel(variableNames)}] = deal('string');
[variableTypes{contains(variableNames,{'date','time'})}] = deal('datetime');
variableTypes(strcmp(variableNames,'squirrel_id')) = {'double'};
% opts = setvartype(opts,opts.VariableNames(contains(opts.VariableNames,{'date','time'})),'datetime');
opts.VariableNames = variableNames;
opts.VariableTypes = variableTypes;

T_AxyDB = readtable(fullfile(storagePath,axyDbFile),opts);
fprintf("total filenames: %i/%i rows\n",sum(strcmp(T_AxyDB.filename,"")==0),size(T_AxyDB,1));

matlabPath = '/Users/matt/Documents/MATLAB/KRSP/Axy';
save(fullfile(matlabPath,'T_AxyDB.mat'),'T_AxyDB','axyDbFile','opts');
%% re-export axy files?
clc

dbFullfiles = string; % compile
for ii = 1:size(T_AxyDB,1)
    dbFullfiles(ii,:) = "";
    if ~ismissing(T_AxyDB.filename(ii))
        csvFiles = dir3(fullfile(storagePath,T_AxyDB.folder{ii}),T_AxyDB.filename{ii});
        if size(csvFiles,1) == 1
            dbFullfiles(ii,:) = string(fullfile(csvFiles.folder,T_AxyDB.filename{ii}));
        else
            fprintf("error %i %s\n",ii,T_AxyDB.filename{ii});
        end
    end
end

axyFullfiles = fullfile(T_AxyFiles.folder,T_AxyFiles.filename);
for ii = 1:numel(dbFullfiles)
%     T_AxyDB.axy_d_fmt(ii) = "";
%     T_AxyDB.axy_t_fmt(ii) = "";
%     T_AxyDB.header_labels(ii) = "";
%     T_AxyDB.axy_fs(ii) = "";

    if ~ismissing(T_AxyDB.filename(ii))
        T_AxyFiles_matchIdx = find(strcmp(axyFullfiles,dbFullfiles(ii)));
        if numel(T_AxyFiles_matchIdx) == 1
            % transfer
            T_AxyDB.axy_d_fmt(ii) = T_AxyFiles.dateFmt(T_AxyFiles_matchIdx);
            T_AxyDB.axy_t_fmt(ii) = T_AxyFiles.timeFmt(T_AxyFiles_matchIdx);
            T_AxyDB.header_labels(ii) = string(strjoin(T_AxyFiles.headerLabels{T_AxyFiles_matchIdx},','));
            T_AxyDB.axy_fs(ii) = string(T_AxyFiles.Fs(T_AxyFiles_matchIdx));
        else
            fprintf("error %i %s\n",ii,T_AxyDB.filename{ii});
        end
    end
end
%% save AxyDB
dateFmt = "dd/MM/yyyy";
timeFmt = "HH:mm:ss";
if isfile(fullfile(storagePath,"._~$"+axyDbFile))
    beep;
    disp("Please close AxyDB before making changes...");
else
    copyfile(fullfile(storagePath,axyDbFile),fullfile(storagePath,axyDbFile+".backup"));
    % can not get this to work dynamically
    T_AxyDB.power_on_date = string(T_AxyDB.power_on_date,dateFmt);
    T_AxyDB.deploy_date = string(T_AxyDB.deploy_date,dateFmt);
    T_AxyDB.removed_date = string(T_AxyDB.removed_date,dateFmt);
    T_AxyDB.power_off_date = string(T_AxyDB.power_off_date,dateFmt);
    T_AxyDB.date_signature = string(T_AxyDB.date_signature,dateFmt);
    T_AxyDB.power_on_time = string(T_AxyDB.power_on_time,timeFmt);
    T_AxyDB.deploy_time = string(T_AxyDB.deploy_time,timeFmt);
    T_AxyDB.removed_time = string(T_AxyDB.removed_time,timeFmt);
    T_AxyDB.power_off_time = string(T_AxyDB.power_off_time,timeFmt);
    T_AxyDB.shake1_start_time = string(T_AxyDB.shake1_start_time,timeFmt);
    T_AxyDB.shake1_stop_time = string(T_AxyDB.shake1_stop_time,timeFmt);
    T_AxyDB.shake2_start_time = string(T_AxyDB.shake2_start_time,timeFmt);
    T_AxyDB.shake2_stop_time = string(T_AxyDB.shake2_stop_time,timeFmt);
    writetable(T_AxyDB,fullfile(storagePath,axyDbFile));
    fprintf("Backup/overwrite of %s successful.\n",axyDbFile);
end

%% find squirrel_id
% krsp_squirrel_file = '/Volumes/GAIDICASSD/KRSP/krsp_squirrel.csv';
% opts = detectImportOptions(krsp_squirrel_file);
% opts = setvartype(opts,opts.VariableNames(contains(opts.VariableNames,{'tag'})),'string');
% krsp_squirrel = readtable(krsp_squirrel_file,opts);
