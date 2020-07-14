if do
    sqkey = readtable('sqkey'); % load meta data table with file locations
    filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp'; % !! change to Temp folder location
    % this loop compiles some data from T for each squirrel
    sq_odba = [];
    sq_awake = [];
    sq_ids = [];
    sq_doys = [];
    sq_dayLength = [];
    iRow = 0;
    squirrelId = 0;
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq}) % make sure .mat file exists for squirrel
            disp(sqkey.filename{iSq}); % we are only looking at this squirrel's data
            load(fullfile(filePath,sqkey.filename{iSq})); % loads T, Tstat
            T = detect_sleepWake(T,2); % adds T.awake column, this function might change
            dtdoys = day(T.datetime,'dayofyear'); % array of days of year
            undoys = unique(dtdoys); % array of unique days of year
            squirrelId = squirrelId + 1; % counts squirrels 
            for iDoy = 1:numel(undoys) % loop through all the days recorded
                theseDoys = find(dtdoys == undoys(iDoy)); % these are the ids in the table for this doy
                if numel(theseDoys) == 1440 % require full day for now (1440 minutes)
                    iRow = iRow + 1; % increment row for these big data arrays
                    sq_ids(iRow) = squirrelId; % log squirrel id, maybe needed later
                    sq_doys(iRow) = undoys(iDoy); % make big array of doys
                    sq_odba(iRow,:) = T.odba(theseDoys); % this should be 1440 data points
                end
            end
        end
    end
    do = false;
end