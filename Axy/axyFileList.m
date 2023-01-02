if do
    tic;
    clc;
    rootDir = '/Volumes/GAIDICASSD/KRSP/KRSP Axy Data';
    files = dir2(rootDir,'-r','*.csv');
    mbSize = 1048576; % MB
    files = files([files.bytes] >= mbSize*50);
    files = struct2table(files);
    
    samplingPeriods = [0.1,1];

    warning ('off','all');
    T_AxyFiles = table;
    for ii = 1:size(files,1)
        T_AxyFiles.id(ii) = ii;
        T_AxyFiles.folder(ii) = files.folder(ii);
        [~,name,ext] = fileparts(files.name(ii));
        T_AxyFiles.filename(ii) = string([name,ext]);
        T_AxyFiles.mb(ii) = round(files.bytes(ii) / mbSize);

        fname = fullfile(T_AxyFiles.folder(ii),T_AxyFiles.filename(ii));
        fid = fopen(fname);
        
        % get fs from first 2 data rows if datetime detected
        for jj = 1:3
            % find delim
            lineText = strrep(fgetl(fid),'"',''); % rep double quotes
            T_AxyFiles.row_header(ii) = {lineText};
            colText = strsplit(lineText,"\t"); % try comma delim
            if numel(colText) == 1
                colText = strsplit(lineText,","); % try tab delim
            end
            if jj == 1 % header
                T_AxyFiles.columns(ii) = numel(colText);
            elseif jj == 2 % data
                T_AxyFiles.row_data(ii) = {lineText};
                dtCols = [];
                for kk = 1:numel(colText)
                    isDt = false;
                    try %#ok<TRYNC>
                        datetime(colText(kk));
                        isDt = true;
                    end
                    if isDt && contains(colText(kk),{'/',':'}) % doubles will register as a dt
                        dtCols = [dtCols kk]; %#ok<AGROW>
                    end
                end
                colTextDiff = colText;
            else
                T_AxyFiles.fs(ii) = NaN;
                for kk = 1:numel(dtCols)
                    diffdt = seconds(datetime(colText(dtCols(kk))) - datetime(colTextDiff(dtCols(kk))));
                    if ismember(diffdt,samplingPeriods)
                        T_AxyFiles.fs(ii) = diffdt;
                        break; % take first match
                    end
                end
            end
        end
        fclose(fid);

        [status,cmdout] = system(sprintf('wc -l "%s"',fname)); % unix
        if status == 0
            cmdParts = strsplit(strtrim(cmdout));
            T_AxyFiles.rows(ii) = str2double(cmdParts(1));
        else
            T_AxyFiles.rows(ii) = NaN;
        end
        if ~isnan(T_AxyFiles.fs(ii))
            T_AxyFiles.days(ii) = floor(T_AxyFiles.rows(ii) / (86400/T_AxyFiles.fs(ii)));
        else
            T_AxyFiles.days(ii) = NaN;
        end

        if any(contains(lower(colText),'temp'))
            T_AxyFiles.has_temp(ii) = true;
        else
            T_AxyFiles.has_temp(ii) = false;
        end
        if any(contains(lower(colText),'odba'))
            T_AxyFiles.has_odba(ii) = true;
        else
            T_AxyFiles.has_odba(ii) = false;
        end

        fprintf("%03d/%i (%1.1f%%) %s - %i days\n",ii,size(files,1),100*ii/size(files,1),T_AxyFiles.filename(ii),T_AxyFiles.days(ii));
    end
    warning ('on','all');
    do = 0;
    toc;
    save("T_AxyFiles",'T_AxyFiles','rootDir');
end
%%
close all;
ff(400,300);

histogram(T_AxyFiles.days);

%% find duplicates
sameList = [];
sameCount = 0;
clc
for ii = 1:size(T_AxyFiles)
    sameIdx = find(T_AxyFiles.rows == T_AxyFiles.rows(ii) & T_AxyFiles.id ~= ii & T_AxyFiles.row_data == T_AxyFiles.row_data(ii));
    for jj = 1:numel(sameIdx)
        if ~ismember(sameIdx(jj),sameList)
            sameCount = sameCount + 1;
            sameList(sameCount) = sameIdx(jj); %#ok<SAGROW>
            fprintf("%s\n%s\n\n",T_AxyFiles.filename(ii),T_AxyFiles.filename(sameIdx(jj)));
        end
    end
end
fprintf("%i/%i entries the same\n",numel(sameList),size(T_AxyFiles,1));

%% find unique headers
unqiueHeaders = unique(T_AxyFiles.row_header);
% for ii = 1:size(T_AxyFiles)
%     if strcmp(T_AxyFiles.row_header,uniqueHeaders
% end