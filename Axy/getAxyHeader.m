function [colMap,colNames,hasHeader,hasError,Fs,startDate,dateFmt,timeFmt,dataLines] = getAxyHeader(fname)
    hasHeader = false;
    hasError = false;
    Fs = NaN;
    dateFmt = "";
    timeFmt = "";
    dateDelim = "/";
    timeDelim = ":";
    startDate = NaT;
    % locale: 'America/Whitehorse'
    
    nFs = 3; % header+2
    fid = fopen(fname);
    dataLines = string;
    jj = 0;
    % gather samples far enough where battery has settled, must include
    % header line (ii=1)
    for ii = 1:2000
        dataLine = fgetl(fid);
        if ii <= nFs || mod(ii,100) == 0
            jj = jj + 1;
            dataLines(jj,:) = strrep(dataLine,'"','');
        end
    end
    % do not close file yet

    colCell = {{'date','Var1'},'time',{'temp','Var5'},...
        {'x','Var2'},{'y','Var3'},{'z','Var4'},'odba',{'battery','alt','Var6'},'nest'};
    colMap = NaN(size(colCell));

    % first look for a labelled header
    if contains(dataLines(1),[colCell{:}],'IgnoreCase',true)
        hasHeader = true;
        headerCols = delimText(dataLines(1));
        for ii = 3:numel(colCell) % skip date/time (below)
            for jj = 1:numel(headerCols)
                if contains(headerCols(jj),[colCell{ii}],'IgnoreCase',true)
                    colMap(ii) = jj;
                    break; % take first match
                end
            end
        end
    end

    % can we detect the rest using data alone?
    dataArr = delimText(dataLines(double(hasHeader)+1:end)); % do not use header
    xlCols = [];
    tempCols = [];
    dateCols = [];
    timeCols = [];
    battCols = [];
    for ii = 1:size(dataArr,2) % iterate columns
        % datetime
        if contains(dataArr(1,ii),dateDelim)
            dateCols = [dateCols ii];
        end
        if contains(dataArr(1,ii),timeDelim)
            timeCols = [timeCols ii];
        end

        dataDouble = str2double(dataArr(:,ii));
        if ~any(isnan(dataDouble)) % this col is a number
            meanCols = mean(dataDouble);
            medCols = median(dataDouble);
            if abs(meanCols) < 3 && abs(medCols) < 2 % x,y,z,odba
                xlCols = [xlCols ii]; %#ok<AGROW>
            elseif medCols > 10
                tempCols = [tempCols ii]; %#ok<AGROW>
            elseif medCols > 3 && medCols < 6
                battCols = [battCols ii]; %#ok<AGROW>
            end
        end
    end

    if ~isempty(dateCols)
        colMap(useCol('date',colCell)) = dateCols(1); % take first
    end
    if ~isempty(timeCols)
        colMap(useCol('time',colCell)) = timeCols(1); % take first
    end
    if isnan(colMap(useCol('x',colCell)))
        if numel(xlCols) > 2
            colMap(useCol('x',colCell)) = xlCols(1);
            colMap(useCol('y',colCell)) = xlCols(2);
            colMap(useCol('z',colCell)) = xlCols(3);
        elseif numel(xlCols) == 1
            % not sure if this condition is in current data
            colMap(useCol('odba',colCell)) = xlCols(1);
        end
    end
    if isnan(colMap(useCol('temp',colCell))) && ~isempty(tempCols)
        colMap(useCol('temp',colCell)) = tempCols(1); % take first
    end
    if isnan(colMap(useCol('battery',colCell))) && ~isempty(battCols)
        colMap(useCol('battery',colCell)) = battCols(1); % take first
    end
    
    % date detector
    yearPos = NaN;
    dayPos = NaN;
    if ~isnan(colMap(useCol('date',colCell)))
        ii = 0;
        while(1)
            dataLine = fgetl(fid);
            if ~isa(dataLine,'double') && mod(ii,1000) == 0 % test periodically
                dataLines(end,:) = strrep(dataLine,'"','');
                dataArr = delimText(dataLines(double(hasHeader)+1:end)); % exclude header
                dateString = dataArr(:,colMap(useCol('date',colCell)));
                try
                    dateParts = split(dateString);
                catch
                    hold on;
                end
                dateNoTime = dateParts(:,1);
                d1 = dateNoTime(end-1);
                d2 = dateNoTime(end);
                if strcmp(d1,d2) == 0 % different
                    d1Parts = strsplit(d1,dateDelim);
                    d2Parts = strsplit(d2,dateDelim);
                    yearPos = find(strlength(d2Parts)==4);
                    possiblePos = 1:3;
                    possiblePos(yearPos) = [];
                    posDiff = [abs(str2double(d1Parts(possiblePos(1))) - str2double(d2Parts(possiblePos(1)))),...
                        abs(str2double(d1Parts(possiblePos(2))) - str2double(d2Parts(possiblePos(2))))];
                    if any(posDiff>0)
                        [~,k] = max(posDiff);
                        dayPos = possiblePos(k);
                        break;
                    end
                end
            end
            ii = ii + 1;
            if isa(dataLine,'double') || ii > 86400*10 % 1day if Fs = 10
                break;
            end
        end
    end
    fclose(fid);
    
    if ~isnan(dayPos) && ~isnan(yearPos)
        dateParts = ["","",""];
        dateParts(dayPos) = "dd";
        dateParts(yearPos) = "yyyy";
        dateParts(strcmp(dateParts,"")) = "MM";
        dateFmt = strjoin(dateParts,"/");
    end
    
    if ~isnan(colMap(useCol('time',colCell)))
        timeString = dataArr(1,colMap(useCol('time',colCell)));
        timeParts = strsplit(timeString,timeDelim);
        if strlength(timeParts(end)) == 2
            timeFmt = "HH:mm:ss";
        else
            timeFmt = "HH:mm:ss.S"; % this should work for all ms after period
        end
    end

    % start error reporting
    if strcmp(dateFmt,"") == 0 && strcmp(timeFmt,"") == 0
        if colMap(useCol('date',colCell)) == colMap(useCol('time',colCell)) % same col
            allDates = datetime(dateString,'inputformat',strjoin([dateFmt,timeFmt]," ")); % join formats
            startDate = allDates(1);
            Fs = 1/seconds(allDates(2)-allDates(1));
        elseif ~isnan(colMap(useCol('time',colCell))) % time is separate
            startDate = datetime(dateString(1),'inputformat',dateFmt);
            timeStrings = dataArr(:,colMap(useCol('time',colCell)));
            allTimes = datetime(timeStrings,'inputformat',timeFmt);
            Fs = 1/seconds(allTimes(2)-allTimes(1));
        end
    end
    if isnan(Fs) || isnat(startDate)
        hasError = true;
        disp("Error: Fs or date issues");
    end

    if isnan(colMap(useCol('odba',colCell)))
        if isnan(colMap(useCol('x',colCell)))
            hasError = true;
            disp("Error: no x/y/z or odba");
        end
    end
    if isnan(colMap(useCol('temp',colCell)))
        hasError = true;
        disp("Error: no temperature");
    end
    
    % compress all possible names
    colNames = colCell;
    for ii = 1:numel(colCell)
        if iscell(colNames{ii})
            colNames(ii) = colNames{ii}(1);
        end
    end
    dataLines = dataLines(1:10); % only return 10
end

function Fs = getFs(data)
    Fs = NaN;
    dateFmt = '';
    try
        dt = datetime(data);
        Fs = 1/seconds(median(diff(dt)));
    catch
        parts = strsplit(data(1),dateDelim);
        if strlength(parts(1)) == 4
            dateFmt = 'yyyy/MM/dd HH:mm:ss';
            dateFmt = dateFmt(1:strlength(data(1))); % only date
            try
                dt = datetime(data,'InputFormat',dateFmt);
                Fs = 1/seconds(median(diff(dt)));
            catch
            end
        end
    end
end