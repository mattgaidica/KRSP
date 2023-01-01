if do
    tic;
    clc;
    rootDir = '/Volumes/GAIDICASSD/KRSP/KRSP Axy Data';
    files = dir2(rootDir,'-r','*.csv');
    mbSize = 1048576; % MB
    files = files([files.bytes] >= mbSize);
    files = struct2table(files);
    samplingPeriods = [1,10];

    warning ('off','all');
    T_AxyFiles = table;
    for ii = 1:size(files,1)
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
            colText = strsplit(lineText,"\t"); % try comma delim
            if numel(colText) == 1
                colText = strsplit(lineText,","); % try tab delim
            end
            if jj == 1 % header
                T_AxyFiles.columns(ii) = numel(colText);
            elseif jj == 2 % data
                dtCols = [];
                for kk = 1:numel(colText)
                    isDt = false;
                    try %#ok<TRYNC>
                        datetime(colText(kk));
                        isDt = true;
                    end
                    if isDt
                        dtCols = [dtCols kk]; %#ok<AGROW>
                    end
                end
                dt = datetime(colText(dtCols));
            else
                T_AxyFiles.fs(ii) = NaN;
                for kk = 1:numel(dtCols)
                    diffdt = seconds(datetime(colText(dtCols(kk))) - dt(kk));
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
end
%%
close all;
ff(400,300);

histogram(T_AxyFiles.days);