% load database
storagePath = '/Volumes/GAIDICASSD/KRSP/Axy Database';
axyDbFile = 'AxyDatabase.xlsx';
opts = detectImportOptions(fullfile(storagePath,axyDbFile));

% override opts
variableNames = opts.VariableNames;
[variableTypes{1:numel(variableNames)}] = deal('string');
[variableTypes{contains(variableNames,{'date','time'})}] = deal('datetime');
[variableTypes{contains(variableNames,{'squirrel_id','axy_fs'})}] = deal('double');
% variableTypes(strcmp(variableNames,'squirrel_id')) = {'double'};
opts.VariableNames = variableNames;
opts.VariableTypes = variableTypes;

T_AxyDB = readtable(fullfile(storagePath,axyDbFile),opts);
fprintf("total filenames: %i/%i rows\n",sum(strcmp(T_AxyDB.filename,"")==0),size(T_AxyDB,1));

%%
savePath = '/Volumes/GAIDICASSD/KRSP/Data_behaviorClass';
hashedAxyFiles = dir3(fullfile(storagePath,'Data'));
iCatch = 0;
iSuccess = 0;
for ii = 1:height(T_AxyDB)
    % find hashed file
    findId = find(contains(hashedAxyFiles.fullfile,T_AxyDB.md5_hash(ii)));
    if ~isempty(findId)
        try
            fprintf("Loading axy data: %s\n",T_AxyDB.md5_hash(ii));
            axy = readtable(hashedAxyFiles.fullfile(findId));
            axy.Properties.VariableNames = strsplit(T_AxyDB.header_labels(ii),",");
            % handle 10Hz data
            if T_AxyDB.axy_fs ~= 1
                % Identify which variables (columns) are of a desired type before
                % downsampling
                validVarTypes = {'double', 'datetime', 'duration'};  % Desired types
                isVarValid = varfun(@(x) any(strcmp(class(x), validVarTypes)), axy, 'OutputFormat', 'uniform');
                % Keep only valid variables
                axy = axy(:,isVarValid);
                axy_timetable = table2timetable(axy);
                resampled_axy = retime(axy_timetable, 'secondly', 'mean');
                axy = timetable2table(resampled_axy);
            end
            % add behavior columns
            disp("filtering temp...");
            axy = april_tempFilter(axy);
            disp("axy prep...");
            axy = april_axyPrep(axy,T_AxyDB.squirrel_id(ii));
            disp("generating statistics...");
            axy = april_summaryStats(axy);
            disp("adding behavior class...");
            axy = april_behavClass2(axy);
            newfilename = sprintf("row%03d_squirrel%i_%s.csv",ii,T_AxyDB.squirrel_id(ii),T_AxyDB.md5_hash(ii));
            fprintf("writing: %s\n",newfilename);
            writetable(axy,fullfile(savePath,newfilename));
            fprintf("Success!\n\n");
            iSuccess = iSuccess + 1;
        catch
            fprintf("ERROR!! %s\n",T_AxyDB.md5_hash(ii));
            iCatch = iCatch + 1;
        end
    end
end
fprintf("Succeeded n = %i\n",iSuccess);
fprintf("Caught n = %i\n",iCatch);