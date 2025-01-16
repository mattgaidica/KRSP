%% Script to Process Squirrel Accelerometer (Axy) Data
% This script processes raw accelerometer data from squirrels and classifies their behaviors.
% It reads from an Excel database containing metadata about each recording and processes
% the corresponding accelerometer files.

% Set up paths and load the database
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

%% Main Processing Loop
% This section processes each accelerometer file in the database:
% 1. Loads the raw data
% 2. Cleans and formats the timestamps
% 3. Removes data outside the deployment period
% 4. Processes the accelerometer signals
% 5. Classifies behaviors
% 6. Saves results and creates summary plots

% Track different types of processing issues
doPlot = true;
skipExisting = false;
savePath = '/Volumes/GAIDICASSD/KRSP/Data_behaviorClass';
hashedAxyFiles = dir3(fullfile(storagePath,'Data'));

% Initialize arrays to track different types of processing issues
iSuccess = 0;
corruptIds = [];      % Files with corrupted data
noDataIds = [];       % Files with no valid data in deployment period
noTempIds = [];       % Files with non-varying temperature
catchIds = [];        % Files that caused errors during processing
successIds = [];      % Successfully processed files

% Process each entry in the database
for ii = 1:height(T_AxyDB)
    % find hashed file
    findId = find(contains(hashedAxyFiles.fullfile,T_AxyDB.md5_hash(ii)));

    if ~isempty(findId)
        % Create a unique filename for the processed data
        newfilename = sprintf("row%03d_squirrel%i_%s",ii,T_AxyDB.squirrel_id(ii),T_AxyDB.md5_hash(ii));
        csvFile = fullfile(savePath,newfilename + '.csv');
        
        % Skip if file already exists and skipExisting is true
        if skipExisting && isfile(csvFile)
            fprintf("Skipping row %i (already exists)\n",ii);
            continue;
        end
        
        try
            fprintf("Loading row %i: %s\n",ii,T_AxyDB.md5_hash(ii));
            % Load the raw accelerometer data
            axy = readtable(hashedAxyFiles.fullfile(findId));

            % Rename columns based on header labels in database
            newHeaderLabels = strsplit(T_AxyDB.header_labels(ii),",");
            axy.Properties.VariableNames(1:numel(newHeaderLabels)) = newHeaderLabels;

            % Check for missing or NaN values in important columns
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

            % Format datetime information consistently
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
            axy.datetime.Format = 'dd-MMM-yyyy HH:mm:ss';

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
            if T_AxyDB.axy_fs(ii) == 10
                % Calculate the subsampling interval based on the sampling frequency
                subsampleInterval = T_AxyDB.axy_fs(ii);

                % Subsample the data
                % Keep the first row and then every 'subsampleInterval'-th row
                axy = axy(1:subsampleInterval:end, :);
            end

            % Process the accelerometer data through several steps:
            disp("filtering temp...");
            axy = april_tempFilter(axy);           % Filter temperature data
            disp("axy prep...");
            axy = april_axyPrep(axy,T_AxyDB.squirrel_id(ii));  % Prepare accelerometer data
            disp("generating statistics...");
            axy = april_summaryStats(axy);         % Calculate summary statistics
            disp("adding behavior class...");
            axy = april_behavClass2(axy);          % Classify behaviors
            fprintf("writing: %s\n",newfilename);

            % Keep only the essential columns in final output
            variables_to_keep = {'datetime','temp','odba','Nest','squirrel','All'};
            variables_to_remove = setdiff(axy.Properties.VariableNames, variables_to_keep);
            axy(:, variables_to_remove) = [];

            % Create visualization if doPlot is true
            if doPlot
                % Subsample data for plotting to improve performance
                numRows = height(axy); % Total number of rows in the data
                maxPoints = 1000;       % Maximum number of points to plot
                interval = max(ceil(numRows / maxPoints), 1);
                idx = (1:interval:numRows)';

                close all;

                % Create a figure with two subplots
                figure('position',[0 0 1200 800]);
                
                % Top subplot: ODBA and Temperature
                subplot(211);
                plot(axy.datetime(idx),axy.odba(idx),'w-');
                set(gca,'ycolor','w');
                xlim([axy.datetime(idx(1)),axy.datetime(idx(end))]);
                xticks(xlim);
                ylabel('ODBA');

                % Add temperature data on secondary y-axis
                yyaxis right;
                plot(axy.datetime(idx),axy.temp(idx),'r-');
                ylabel('Temp (C)');
                set(gca,'ycolor','r');
                hold on;
                
                % Add nest occupancy indicator
                nestData = double(contains(axy.Nest(idx),"Out"));
                nestData(nestData == 1) = NaN;
                plot(axy.datetime(idx),nestData + mean(ylim),'g-','LineWidth',3);
                legend({"ODBA","Temp","In Nest"});
                title('Axy + Temp')

                % Bottom subplot: Behavior Classification
                subplot(212);
                % Convert behavior classes to categorical and set specific order
                allCategorical = categorical(axy.All);
                desiredOrder = {'NestNotMove', 'NestMove', 'NotMoving', 'Feed', 'Forage', 'Travel'};
                allCategorical = reordercats(allCategorical, desiredOrder);

                % Plot behavior classes over time
                plot(axy.datetime(idx),allCategorical(idx));
                xlim([axy.datetime(idx(1)),axy.datetime(idx(end))]);
                xticks(xlim);

                % Add nest occupancy indicator
                yyaxis right;
                plot(axy.datetime(idx),nestData + mean(ylim),'g-','LineWidth',3);
                set(gca,'ycolor','w');
                yticks([]);
                title('Behavior Class');
                legend({"Class","In Nest"});

                % Save the figure as a JPG
                exportgraphics(gcf,fullfile(savePath,newfilename + '.jpg'));
            end

            % Save processed data to CSV
            writetable(axy,csvFile);
            fprintf("Success!\n\n");
            successIds = [successIds;ii];
        catch
            fprintf("!! Caught %s\n",T_AxyDB.md5_hash(ii));
            catchIds = [catchIds;ii];
        end
    end
end

%% Final Summary
% Display counts of different processing outcomes
fprintf("Caught n = %i\n",numel(catchIds));
fprintf("Corrupt n = %i\n",numel(corruptIds));
fprintf("No data n = %i\n",numel(noDataIds));
fprintf("No temp n = %i\n",numel(noTempIds));