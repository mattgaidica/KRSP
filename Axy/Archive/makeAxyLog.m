storagePath = '/Volumes/GAIDICASSD/KRSP';
axyPath = '/Users/matt/Documents/MATLAB/KRSP/Axy';
axyDbFile = '/Volumes/GAIDICASSD/KRSP/KRSP Axy Data/AxyDatabase.xlsx';
opts = detectImportOptions(axyDbFile);

opts = setvartype(opts,opts.VariableNames(contains(opts.VariableNames,{'shake','time'})),'datetime');
T_AxyDB = readtable(axyDbFile,opts);
fprintf("total filenames: %i/%i rows\n",sum(strcmp(T_AxyDB.filename,"")==0),size(T_AxyDB,1));
save(fullfile(axyPath,'T_AxyDB.mat'),'T_AxyDB','storagePath','axyDbFile','opts');

%%
rowIds = find(strcmp(T_AxyDB.logfile,'SQRaxy_key_mid_den.xlsx'));
csvFiles = dir2(fullfile(storagePath,T_AxyDB.root_dir{rowIds(1)}),'-r','*.csv');
[~,filenames] = fileparts({csvFiles.name});

tempFiles = string;
fileDeployments = string;
fileMiddens = string;
jj = 0;
for ii = 1:numel(filenames)
    if strcmp(filenames{ii}(1),'.') || csvFiles(ii).bytes < 100000 % hidden
        continue;
    end
    jj = jj + 1;
    tempFiles(jj,1) = string(filenames{ii})+'.csv';
    csvParts = split(csvFiles(ii).name,filesep);
    fileDeployments(jj,1) = strrep(string(csvParts{1})," - "," ");
    filenameParts = split(tempFiles(jj,1),"_");
    a = filenameParts{1};
    a = strrep(a,".","");
    a = strrep(a,"-","n");
    a = strrep(a,"j","");
    fileMiddens(jj,1) = string(a);
end
filenames = tempFiles;

clc
matchFiles = string;
for ii = 1:numel(rowIds)
    a = T_AxyDB.midden{rowIds(ii)};
    a = strrep(a,".","");
    a = strrep(a,"-","n");
    matchIdx = find(strcmp(fileMiddens,string(a)) & strcmp(fileDeployments,T_AxyDB.deployment{rowIds(ii)}));
    if ~isempty(matchIdx)
        if numel(matchIdx) == 1
            matchFiles(ii,1) = filenames{matchIdx};
        else
            fprintf("multiple %i row %",ii,rowIds(ii)+1);
        end
    else
        matchFiles(ii,1) = "";
        fprintf("%i row %i not found\n",ii,rowIds(ii)+1);
    end
end
%% needs home
clc
for ii = 1:numel(matchFiles)
    if ~isempty(matchFiles{ii})
        matchIdx = find(strcmp(matchFiles,matchFiles{ii}));
        if numel(matchIdx) > 1
            for jj = 1:numel(matchIdx)
                fprintf("duplicate %s at %i (row %i) \n",matchFiles{matchIdx(jj)},matchIdx(jj),rowIds(matchIdx(jj))+1);
            end
            fprintf("\n");
        end
    end
end

for ii = 1:numel(filenames)
    matchIdx = find(strcmp(matchFiles,filenames{ii}));
    if isempty(matchIdx)
        csvId = find(contains({csvFiles.name},filenames{ii}));
        fprintf("%s needs home\n",csvFiles(csvId).name);
    end
end

%%
rowIds = find(strcmp(T_AxyDB.logfile,'2019 Playback Experiment-Data_16Oct.xlsx'));

useMiddens = "";
for ii = 1:numel(modexport(:,2))
    if ~isempty(modexport{ii,2})
        useMiddens(ii,1) = strrep(modexport{ii,2},".","");
    else
        useMiddens(ii,1) = "";
    end
end


matchFiles = string;
for ii = 1:numel(rowIds)
   thisMidden = strrep(T_AxyDB.midden{rowIds(ii)},".","");
   matchIdx = find(strcmp(modexport(:,1),T_AxyDB.grid{rowIds(ii)}) & strcmp(useMiddens,thisMidden) &...
       transpose([modexport{:,4}]==T_AxyDB.squirrel_id(rowIds(ii))) & transpose([modexport{:,6}]==str2double(T_AxyDB.axy{rowIds(ii)})) &...
       T_AxyDB.power_on_date(rowIds(ii)) == a);
   if numel(matchIdx) == 1 && ~isempty(modexport{matchIdx,5})
        matchFiles(ii,1) = modexport{matchIdx,5};
   else
        matchFiles(ii,1) = "";
   end
end
%%
% % rowIds = find(strcmp(T_AxyDB.logfile,'2019 Playback Experiment-Data_16Oct.xlsx'));
% % csvFiles = dir2(fullfile(storagePath,T_AxyDB.root_dir{rowIds(1)}),'-r','*.csv');
% % [~,filenames] = fileparts({csvFiles.name});
% % 
% % jj = 0;
% % tempFiles = string;
% % axys = string;
% % squirrelIds = [];
% % for ii = 1:numel(filenames)
% %     csvParts = strsplit(csvFiles(ii).name,filesep);
% %     if strcmp(filenames{ii}(1),'.') || csvFiles(ii).bytes < 100000 % hidden/kB
% %         continue;
% %     end
% %     jj = jj + 1;
% %     tempFiles(jj,1) = string(filenames{ii})+'.csv';
% %     nameParts = strsplit(filenames{ii},'_');
% %     axys(jj,:) = nameParts{1}(4:end);
% % end
% % filenames = tempFiles;
% % % % % % %% % get these from old file
% % % % % % squirreIds = [];
% % % % % % for ii = 1:numel(rowIds)
% % % % % %     tagIds = find(strcmp(T_AxyDB.tag_left{rowIds(ii)},modexport(:,2)) & strcmp(T_AxyDB.tag_right{rowIds(ii)},modexport(:,3)));
% % % % % %     if ~isempty(tagIds)
% % % % % %         squirreIds(ii,1) = modexport{tagIds,1};
% % % % % %     else
% % % % % %         squirreIds(ii,1) = NaN;
% % % % % %     end
% % % % % % end

clc
matchFiles = string;
for ii = 1:numel(rowIds)
    matchIdx = find(strcmp(axys,T_AxyDB.axy{rowIds(ii)}) & contains(filenames,string(T_AxyDB.squirrel_id(rowIds(ii)))));
    if numel(matchIdx) > 1
        fprintf("%i !! MULTIPLE\n",rowIds(ii)+1);
    end
    if numel(matchIdx) == 1
        matchFiles(ii,1) = sprintf("%s",filenames{matchIdx(1)});
    end
    if numel(matchIdx) > 1
        for jj = 1:numel(matchIdx)
           fprintf("%s\n",filenames{matchIdx(jj)});
        end
        fprintf("\n");
        matchFiles(ii,1) = "";
    end
    if isempty(matchIdx)
        matchFiles(ii,1) = "";
    end
end
%%
clc
for ii = 1:numel(matchFiles)
    if ~isempty(matchFiles{ii})
        matchIdx = find(strcmp(matchFiles,matchFiles{ii}));
        if numel(matchIdx) > 1
            for jj = 1:numel(matchIdx)
                fprintf("duplicate %s at %i (row %i) \n",matchFiles{matchIdx(jj)},matchIdx(jj),rowIds(matchIdx(jj))+1);
            end
            fprintf("\n");
        end
    end
end

for ii = 1:numel(filenames)
    matchIdx = find(strcmp(matchFiles,filenames{ii}));
    if isempty(matchIdx)
        fprintf("%s needs home\n",filenames{ii});
    end
end

%%
rowIds = find(strcmp(T_AxyDB.logfile,'AxyLog2017.xlsx'));
csvFiles = dir2(fullfile(storagePath,T_AxyDB.root_dir{rowIds(1)}),'-r','*.csv');
[~,filenames] = fileparts({csvFiles.name});
jj = 0;

tempFiles = string;
for ii = 1:numel(filenames)
    if strcmp(filenames{ii}(1),'.') % hidden
        continue;
    end
    jj = jj + 1;
    tempFiles(jj,1) = string(filenames{ii})+'.csv';
end
filenames = tempFiles;

clc
matchFiles = string;
for ii = 1:numel(rowIds)
    matchIdx = find(contains(filenames,T_AxyDB.axy{rowIds(ii)}) & contains(filenames,T_AxyDB.grid{rowIds(ii)}) & ...
    contains(filenames,string(T_AxyDB.squirrel_id(rowIds(ii)))));
    if numel(matchIdx) > 1
        fprintf("%i !! MULTIPLE\n",rowIds(ii)+1);
    end
    if ~isempty(matchIdx)
        matchFiles(ii,1) = sprintf("%s",filenames{matchIdx(1)});
    else
        disp("empty");
        matchFiles(ii,1) = "";
    end
end

%%
clc
for ii = 1:numel(matchFiles)
    if ~isempty(matchFiles{ii})
        matchIdx = find(strcmp(matchFiles,matchFiles{ii}));
        if numel(matchIdx) > 1
            for jj = 1:numel(matchIdx)
                fprintf("duplicate %s \n",matchFiles{matchIdx(jj)});
            end
            fprintf("\n");
        end
    end
end

for ii = 1:numel(filenames)
    matchIdx = find(strcmp(matchFiles,filenames{ii}));
    if isempty(matchIdx)
        fprintf("%s needs home\n",filenames{ii});
    end
end
%% what a mess
rowIds = find(strcmp(T_AxyDB.logfile,'AxyLogSquirrel-2018.xlsx'));
csvFiles = dir2(fullfile(storagePath,T_AxyDB.root_dir{rowIds(1)}),'-r','*.csv');
[~,filenames] = fileparts({csvFiles.name});

jj = 0;
searchFilenames = {};
tempFiles = {};
for ii = 1:numel(filenames)
    if strcmp(filenames{ii}(1),'.') % hidden
        continue;
    end
    jj = jj + 1;
    tempFiles{jj} = filenames{ii};
    parts = strsplit(filenames{ii},'_');
    searchFilenames{jj} = parts{1};
end
filenames = tempFiles;

clc
matchFiles = string;
for ii = 1:numel(rowIds)
    matchIdx = find((strcmp(searchFilenames,T_AxyDB.axy{rowIds(ii)}) & contains(filenames,T_AxyDB.grid{rowIds(ii)})) | ...
        contains(filenames,string(T_AxyDB.squirrel_id(rowIds(ii)))));
    if ~isempty(matchIdx)
        if numel(matchIdx) > 1
            fprintf("%i !! MULTIPLE\n",rowIds(ii)+1);
            % more specific
%             matchIdx = find((strcmp(searchFilenames,T_AxyDB.axy{rowIds(ii)}) & contains(filenames,T_AxyDB.grid{rowIds(ii)})));
%             break;
        end
        for jj = 1:numel(matchIdx)
            thisFilename = sprintf("%s.csv",filenames{matchIdx(jj)});
            axyFileId = find(strcmp(T_AxyFiles.filename,thisFilename));
            fprintf("%i (%i) match for %s.csv\n",rowIds(ii)+1,ii,filenames{matchIdx(jj)});
        end
        if numel(matchIdx) == 1
            matchFiles(ii,1) = sprintf("%s.csv",filenames{matchIdx(1)});
        else
            matchFiles(ii,1) = "";
        end
        fprintf("\n");
    else
        fprintf("!! %i (%i) no match\n",rowIds(ii)+1,ii);
        matchFiles(ii,1) = "";
    end
end

clc
foundIds = [];
for ii = 1:numel(matchFiles)
    foundFiles = find(strcmp(filenames,matchFiles{ii}(1:end-4)));
    if numel(foundFiles)>0 && ~isempty(matchFiles{ii})
        foundIds = [foundIds foundFiles];
    end

%     if ~isempty(matchFiles{ii})
%         matchIdx = find(strcmp(matchFiles,matchFiles{ii}));
%         if numel(matchIdx) > 1
%             for jj = 1:numel(matchIdx)
%                 fprintf("%s \n",matchFiles{matchIdx(jj)});
%             end
%             fprintf("\n");
%         end
%     end
end
%%
rowIds = find(strcmp(T_AxyDB.logfile,'AxyLogSquirrel-2019.xlsx'));
csvFiles = dir2(fullfile(storagePath,T_AxyDB.root_dir{rowIds(1)}),'-r','*.csv');
[~,filenames] = fileparts({csvFiles.name});
filenames(contains(filenames,'.')) = [];
searchFilenames = filenames;
for ii = 1:numel(filenames)
    parts = strsplit(filenames{ii},'_');
    searchFilenames{ii} = parts{1};
end
clc
matchFiles = string;
for ii = 1:numel(rowIds)
    matchIdx = find(strcmp(searchFilenames,T_AxyDB.axy{rowIds(ii)}));
    if ~isempty(matchIdx)
        if numel(matchIdx) > 1
            fprintf("%i !! MULTIPLE\n",ii);
            break;
        end
        thisFilename = sprintf("%s.csv",filenames{matchIdx(1)});
        axyFileId = find(strcmp(T_AxyFiles.filename,thisFilename));
        fprintf("%i match for %s.csv\n",ii,filenames{matchIdx(1)});
        matchFiles(ii,1) = sprintf("%s.csv",filenames{matchIdx(1)});
    else
        fprintf("!! %i no match\n",ii);
        matchFiles(ii,1) = "";
    end
end

%%
rowIds = find(strcmp(T_AxyDB.logfile,'JO Axy Log 2016.xlsx'));
csvFiles = dir2(fullfile(storagePath,T_AxyDB.root_dir{rowIds(1)}),'-r','*.csv');
[~,filenames] = fileparts({csvFiles.name});
clc
matchFiles = string;
matchIds = [];
for ii = 1:numel(rowIds)
    treatment = T_AxyDB.treatment{ii}; % GC
    if strcmp(T_AxyDB.treatment{ii},'CONTROL')
        treatment = 'ctrl';
    elseif strcmp(T_AxyDB.treatment{ii},'NONE')
        treatment = 'none';
    end
    breeding = lower(T_AxyDB.breeding_status{ii}); % lactation
    if strcmp(T_AxyDB.breeding_status{ii},'PREGNANCY')
        breeding = 'preg';
    elseif strcmp(T_AxyDB.breeding_status{ii},'NON BREEDING')
        breeding = 'nb';
    end
    matchIdx = find(contains(filenames,string(T_AxyDB.squirrel_id(rowIds(ii)))) &...
        contains(filenames,T_AxyDB.midden{rowIds(ii)}) &...
        contains(filenames,treatment) &...
        contains(filenames,breeding));
    if ~isempty(matchIdx)
        if numel(matchIdx) > 1
            fprintf("%i !! MULTIPLE\n",ii);
            break;
        end
        thisFilename = sprintf("%s.csv",filenames{matchIdx(1)});
        axyFileId = find(strcmp(T_AxyFiles.filename,thisFilename));
        startDate = T_AxyFiles.startDate(axyFileId);
%         fprintf("%i match for %s.csv | %s (%i) - %s\n",ii,filenames{matchIdx(jj)},string(startDate),axyFileId,T_AxyDB.deployment{ii});
        matchFiles(ii,1) = sprintf("%s.csv",filenames{matchIdx(1)});
        matchIds = [matchIds matchIdx];
    else
        fprintf("!! %i no match\n",ii);
        matchFiles(ii,1) = "";
    end
%     fprintf("\n");
end
noMatch = ~ismember(matchIds,1:numel(filenames));
disp(noMatch);
%%
rowIds = find(strcmp(T_AxyDB.logfile,'AxyLog2016.xlsx'));
csvFiles = dir2(fullfile(storagePath,T_AxyDB.root_dir{rowIds(1)}),'-r','*.csv');
[~,filenames] = fileparts({csvFiles.name});
filenames(contains(filenames,'.')) = [];
clc
matchFiles = string;
for ii = 1:numel(rowIds)
    matchIdx = find(strcmp(filenames,[T_AxyDB.axy{rowIds(ii)},'_1']));
    if ~isempty(matchIdx)
        if numel(matchIdx) > 1
            fprintf("%i !! MULTIPLE\n",ii);
            break;
        end
        thisFilename = sprintf("%s.csv",filenames{matchIdx(1)});
        axyFileId = find(strcmp(T_AxyFiles.filename,thisFilename));
        fprintf("%i match for %s.csv\n",ii,filenames{matchIdx(jj)});
        matchFiles(ii,1) = sprintf("%s.csv",filenames{matchIdx(1)});
    else
        fprintf("!! %i no match\n",ii);
        matchFiles(ii,1) = "";
    end
end