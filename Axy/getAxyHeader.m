function [colMap,colNames,hasHeader,hasError,Fs,startDate,dateFmt,dataLines] = getAxyHeader(fname)
    hasHeader = false;
    hasError = false;
    Fs = NaN;
    dateFmt = '';
    startDate = NaT;
    
    fid = fopen(fname);
    dataLines = string;
    jj = 0;
    for ii = 1:2000
        dataLine = fgetl(fid);
        if mod(ii,200) == 0
            jj = jj + 1;
            dataLines(jj,:) = strrep(dataLine,'"','');
        end
    end
    fclose(fid);

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
        if contains(dataArr(1,ii),'/')
            dateCols = [dateCols ii];
        end
        if contains(dataArr(1,ii),':')
            timeCols = [timeCols ii];
        end

        dataDouble = str2double(dataArr(:,ii));
        if ~any(isnan(dataDouble)) % this col is a number
            meanCols = mean(dataDouble);
            medCols = median(dataDouble);
            if abs(meanCols) < 3 && abs(medCols) < 2 % x,y,z,odba
                xlCols = [xlCols ii]; %#ok<AGROW>
            elseif meanCols > 10
                tempCols = [tempCols ii]; %#ok<AGROW>
            elseif meanCols > 3 && meanCols < 6
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
    if colMap(useCol('date',colCell)) == colMap(useCol('time',colCell)) % same col
        [Fs,dateFmt] = getFs(dataArr(:,colMap(useCol('date',colCell))));
    elseif ~isnan(colMap(useCol('time',colCell))) % time is separate
        Fs = getFs(dataArr(:,colMap(useCol('time',colCell)))); % Fs from time
        [~,dateFmt] = getFs(dataArr(:,colMap(useCol('date',colCell)))); % fmt from date
    end
    % combine to get startDate
    if ~isnan(colMap(useCol('date',colCell)))
        dateString = dataArr(1,colMap(useCol('date',colCell)));
        dateParts = strsplit(dateString);
        if isempty(dateFmt)
            startDate = datetime(dateParts(1));
        else
            fmtParts = strsplit(dateFmt); % only return date (not time)
            startDate = datetime(dateParts(1),'inputformat',fmtParts{1});
        end
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

    % compress all possible names
    colNames = colCell;
    for ii = 1:numel(colCell)
        if iscell(colNames{ii})
            colNames(ii) = colNames{ii}(1);
        end
    end

    % can we form a datetime?
    if isnan(colMap(useCol('datetime',colCell)))
        if isnan(colMap(useCol('date',colCell))) || isnan(colMap(useCol('time',colCell)))
            hasError = true;
            disp("Error: unable to make datetime");
        end
    end
    if isnan(Fs)
        hasError = true;
        disp("Error: Fs was not determined");
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
end

function [Fs,dateFmt] = getFs(data)
    Fs = NaN;
    dateFmt = '';
    try
        dt = datetime(data);
        Fs = 1/seconds(median(diff(dt)));
    catch
        parts = strsplit(data(1),'/');
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

function dataArr = delimText(lineText)
    for ii = 1:size(lineText,1)
        colText = strsplit(lineText(ii),"\t"); % try tab delim
        if numel(colText) == 1 % try again
            colText = strsplit(lineText(ii),","); % try comma delim
        end
        dashIdx = find(count(colText,'-') > 1); % must be a date
        colText(dashIdx) = strrep(colText(dashIdx),'-','/');
        dataArr(ii,:) = colText; %#ok<AGROW>
    end
end