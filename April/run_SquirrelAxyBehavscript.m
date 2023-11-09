% upload: https://drive.google.com/drive/folders/1z37UOIedbWeJea2c4oOftKt6ghXA8eOV
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
% check: c92f49590fdfd9a03cc5df6232a27268, 5a8d427fb949f28ac5cfa67d1bceb83b
doPlot = true;
skipExisting = true;
savePath = '/Volumes/GAIDICASSD/KRSP/Data_behaviorClass';
hashedAxyFiles = dir3(fullfile(storagePath,'Data'));
iSuccess = 0;
corruptIds = [];
noDataIds = [];
noTempIds = [];
catchIds = [];
successIds = [];
for ii = 676:height(T_AxyDB)
    % find hashed file
    findId = find(contains(hashedAxyFiles.fullfile,T_AxyDB.md5_hash(ii)));

    if ~isempty(findId)
        newfilename = sprintf("row%03d_squirrel%i_%s",ii,T_AxyDB.squirrel_id(ii),T_AxyDB.md5_hash(ii));
        csvFile = fullfile(savePath,newfilename + '.csv');
        if skipExisting && isfile(csvFile)
            fprintf("Skipping row %i (already exists)\n",ii);
            continue;
        end
        try
            fprintf("Loading row %i: %s\n",ii,T_AxyDB.md5_hash(ii));
            % not using axy_d/t_format here, maybe I should?
            axy = readtable(hashedAxyFiles.fullfile(findId)); % ,'Delimiter',delimiter

            newHeaderLabels = strsplit(T_AxyDB.header_labels(ii),",");
            axy.Properties.VariableNames(1:numel(newHeaderLabels)) = newHeaderLabels;

            containsNaNorMissing = false;
            variables_to_check = {'datetime','date','time','temp','odba','x','y','z'};
            for col = 1:width(axy)
                if contains(axy.Properties.VariableNames{col},variables_to_check)
                    columnData = axy{:,col};
                    if isnumeric(columnData) && any(isnan(columnData))
                        containsNaNorMissing = true;
                    end
                end
            end
            if containsNaNorMissing
                disp("Data corrupt, skipping");
                corruptIds = [corruptIds;ii];
                continue;
            end

            % handle case where date/time is string and combine them
            t_fmt = T_AxyDB.axy_t_fmt(ii);
            d_fmt = T_AxyDB.axy_d_fmt(ii);
            if ~any(contains(axy.Properties.VariableNames,"datetime"))
                if iscell(axy.date)
                    if sum(axy.date{1}=='-') == 2
                        d_fmt = strrep(d_fmt,'/','-'); % replace
                    end
                    axy.date = datetime(axy.date,'Format',d_fmt);
                end
                if iscell(axy.time)
                    axy.time = datetime(axy.time,'Format',t_fmt);
                end
                % combine date, time
                axy.datetime = axy.date + axy.time;
            end
            % !! fixing similar issue to above but requires datetime
            datetimeIncorrect = false;
            if iscell(axy.datetime)
                try
                    if sum(axy.datetime{1}=='-') == 2
                        d_fmt = strrep(T_AxyDB.axy_d_fmt(ii),'/','-'); % replace
                    end
                    axy.datetime = datetime(axy.datetime,'Format',d_fmt+" "+t_fmt);
                catch
                    datetimeIncorrect = true;
                end
            end
            if datetimeIncorrect
                disp("failed to convert datetime from string, skipping");
                continue;
            end

            % remove nose/tail of data
            % handle the way date/time is read
            disp("removing non-deployment periods...")
            noseDatetime = datetime(year(T_AxyDB.deploy_date(ii)),month(T_AxyDB.deploy_date(ii)),day(T_AxyDB.deploy_date(ii)),...
                hour(T_AxyDB.deploy_time(ii)),minute(T_AxyDB.deploy_time(ii)),second(T_AxyDB.deploy_time(ii)));
            tailDatetime = datetime(year(T_AxyDB.removed_date(ii)),month(T_AxyDB.removed_date(ii)),day(T_AxyDB.removed_date(ii)),...
                hour(T_AxyDB.removed_time(ii)),minute(T_AxyDB.removed_time(ii)),second(T_AxyDB.removed_time(ii)));
            % check for valid times, keep going if it can't trim (need to
            % review plot to find valid deployment interval)
            if ~isnat(noseDatetime) && ~isnat(tailDatetime)
                axy(axy.datetime < noseDatetime | axy.datetime > tailDatetime,:) = []; % remove
            end
            if isempty(axy)
                disp("No data in range, skipping");
                noDataIds = [noDataIds;ii];
                continue;
            end
            if std(axy.temp) == 0
                disp("Non-varying temperature, skipping");
                noTempIds = [noTempIds;ii];
                continue;
            end

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
            fprintf("writing: %s\n",newfilename);

            % List of variables to keep
            variables_to_keep = {'datetime','temp','odba','Nest','squirrel','All'};
            % Find variables in the table that are not in the list to keep
            variables_to_remove = setdiff(axy.Properties.VariableNames, variables_to_keep);
            % Remove unwanted variables
            axy(:, variables_to_remove) = [];

            % create summary figure
            if doPlot
                numRows = height(axy); % Total number of rows in the data
                maxPoints = 1000;       % Maximum number of points to plot
                interval = max(ceil(numRows / maxPoints), 1);
                idx = (1:interval:numRows)';
                
                close all;
                
                figure('position',[0 0 1200 800]);
                subplot(211);
                plot(axy.datetime(idx),axy.odba(idx),'w-');
                set(gca,'ycolor','w');
                xlim([axy.datetime(idx(1)),axy.datetime(idx(end))]);
                xticks(xlim);
                ylabel('ODBA');
                
                yyaxis right;
                plot(axy.datetime(idx),axy.temp(idx),'r-');
                ylabel('Temp (C)');
                set(gca,'ycolor','r');
                hold on;
                nestData = double(contains(axy.Nest(idx),"Out"));
                nestData(nestData == 1) = NaN;
                plot(axy.datetime(idx),nestData + mean(ylim),'g-','LineWidth',3);
                legend({"ODBA","Temp","In Nest"});
                title('Axy + Temp')
                
                subplot(212);
                plot(axy.datetime(idx),categorical(axy.All(idx)));
                xlim([axy.datetime(idx(1)),axy.datetime(idx(end))]);
                xticks(xlim);
                
                yyaxis right;
                plot(axy.datetime(idx),nestData + mean(ylim),'g-','LineWidth',3);
                set(gca,'ycolor','w');
                yticks([]);
                title('Behavior Class');
                legend({"Class","In Nest"});

                exportgraphics(gcf,fullfile(savePath,newfilename + '.jpg'));
            end

            % write CSV file
            writetable(axy,csvFile);
            fprintf("Success!\n\n");
            successIds = [successIds;ii];
        catch
            fprintf("!! Caught %s\n",T_AxyDB.md5_hash(ii));
            catchIds = [catchIds;ii];
        end
    end
end

%%
fprintf("Caught n = %i\n",numel(catchIds));
fprintf("Corrupt n = %i\n",numel(corruptIds));
fprintf("No data n = %i\n",numel(noDataIds));
fprintf("No temp n = %i\n",numel(noTempIds));