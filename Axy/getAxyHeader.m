function results = getAxyHeader(fname)
hasHeader = false;
hasError = false;
Fs = NaN;
dateFmt = "";
timeFmt = "";
dateDelim = "/";
timeDelim = ":";
startDate = NaT;
nRows = NaN;
nDays = NaN;
headerCols = {};
md5 = "";
% locale: 'America/Whitehorse'

nFs = 5; % header+2
fid = fopen(fname);
dataLines = string;
jj = 0;
% gather samples far enough where battery has settled, must include
% header line (ii=1)
for ii = 1:2000
    dataLine = fgetl(fid);
    if ~isa(dataLine,'double') && (ii <= nFs || mod(ii,100) == 0)
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
        dateCols = [dateCols ii]; %#ok<AGROW> 
    end
    if contains(dataArr(1,ii),timeDelim)
        timeCols = [timeCols ii]; %#ok<AGROW> 
    end

    dataDouble = str2double(dataArr(:,ii));
    if ~any(isnan(dataDouble)) % this col is a number
        meanCols = mean(dataDouble);
        medCols = median(dataDouble);
        if abs(meanCols) < 3 && abs(medCols) < 2 % x,y,z,odba
            xlCols = [xlCols ii]; %#ok<AGROW>
        elseif medCols > 10
            tempCols = [tempCols ii]; %#ok<AGROW>
        else
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
        if mod(ii,7200) == 0 % test periodically
            if ~isa(dataLine,'double') % in case -1 from reading file above
                dataLines(end,:) = strrep(dataLine,'"','');
                dataArr = delimText(dataLines(double(hasHeader)+1:end)); % exclude header
            end
            dateString = dataArr(:,colMap(useCol('date',colCell)));
            dateParts = split(dateString);
            dateNoTime = dateParts(:,1);
            d1 = dateNoTime(end-1);
            d2 = dateNoTime(end);
            d1Parts = double(strsplit(d1,dateDelim));
            d2Parts = double(strsplit(d2,dateDelim));
            yearPos = find(d1Parts>2000);
            if isempty(yearPos) % year must not exist/two digits
                break;
            end
            possiblePos = 1:3;
            possiblePos(yearPos) = [];
            % detect day by > 12
            k = find(d1Parts(possiblePos)>12 | d2Parts(possiblePos)>12);
            if ~isempty(k)
                dayPos = possiblePos(k);
                break;
            end
            % or detect day by changing values
            if strcmp(d1,d2) == 0 % different
                posDiff = abs(d1Parts(possiblePos)-d2Parts(possiblePos));
                if any(posDiff>0)
                    [~,k] = max(posDiff);
                    dayPos = possiblePos(k);
                    break;
                end
            end
        end
        % read line after above has executed once
        dataLine = fgetl(fid);
        if isa(dataLine,'double') || ii > 86400*10 % 1day if Fs = 10
            break;
        end
        ii = ii + 1;
    end
end
fclose(fid);

if isnan(dayPos) && ~isnan(yearPos)
    if yearPos == 1
        dayPos = 3;
    elseif yearPos == 3
        dayPos = 1;
    end
end

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

% % % small files?
% % if strcmp(dateFmt,"")
% %     if colMap(useCol('date',colCell)) == 1
% %         dateFmt = "dd/MM/yyyy";
% %     end
% % end

% start error reporting
if strcmp(dateFmt,"") == 0 && strcmp(timeFmt,"") == 0
    dateString = dataArr(:,colMap(useCol('date',colCell))); % make sure this exists
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

[status,cmdout] = system(sprintf('wc -l "%s"',fname)); % unix
if status == 0
    cmdParts = strsplit(strtrim(cmdout));
    nRows = str2double(cmdParts(1));
end

[status,cmdout] = system(sprintf('md5 "%s"',fname)); % unix
if status == 0
    cmdParts = strsplit(strtrim(cmdout));
    md5 = string(cmdParts(end));
end

if ~isnan(nRows) && ~isnan(Fs)
    nDays = floor(nRows / (86400*Fs));
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
% make header labels
headerLabels = cell(1,size(dataArr,2));
for ii = 1:numel(colMap)
    if ~isnan(colMap(ii))
        if ~isempty(headerLabels{colMap(ii)}) % datetime condition
            headerLabels{colMap(ii)} = [headerLabels{colMap(ii)},colNames{ii}];
        else
            headerLabels{colMap(ii)} = colNames{ii};
        end
    end
end

% fill missing
defaultLabel  = "var";
defaultCount = 0;
for ii = 1:numel(headerLabels)
    if isempty(headerLabels{ii})
        defaultCount = defaultCount + 1;
        headerLabels{ii} = sprintf('%s%i',defaultLabel,defaultCount);
    end
end

% make results struct, !!could replace inline
results = struct;
results.colMap = colMap;
results.colNames = colNames;
results.headerColumns = numel(headerCols);
results.headerLabels = headerLabels;
results.hasHeader = hasHeader;
results.hasError = hasError;
results.Fs = Fs;
results.startDate = startDate;
results.dateFmt = dateFmt;
results.timeFmt = timeFmt;
results.dataLines = dataLines(1:nFs);
results.days = nDays;
results.rows = nRows;
results.md5 = md5;
