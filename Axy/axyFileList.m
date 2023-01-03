if do
    tic;
    clc;
    rootDir = '/Volumes/GAIDICASSD/KRSP/KRSP Axy Data';
    files = dir2(rootDir,'-r','*.csv');
    mbSize = 1048576; % MB
    files = files([files.bytes] >= mbSize*30);
    files = struct2table(files);

    warning ('off','all');
    T_AxyFiles = table;
    for ii = 1:size(files,1)
        T_AxyFiles.id(ii) = ii;
        T_AxyFiles.folder(ii) = string(files.folder(ii));
        [~,name,ext] = fileparts(files.name(ii));
        T_AxyFiles.filename(ii) = string([name,ext]);
        T_AxyFiles.mb(ii) = round(files.bytes(ii) / mbSize);

        fname = fullfile(T_AxyFiles.folder(ii),T_AxyFiles.filename(ii));
        
        [colMap,colNames,hasHeader,hasError,Fs,startDate,dateFmt,timeFmt,dataLines] = getAxyHeader(fname);
        T_AxyFiles.Fs(ii) = Fs;
        T_AxyFiles.colNames(ii) = {colNames};
        T_AxyFiles.colMap(ii) = {colMap};
        T_AxyFiles.hasError(ii) = hasError;
        T_AxyFiles.hasHeader(ii) = hasHeader;
        T_AxyFiles.startDate(ii) = startDate;
        T_AxyFiles.dateFmt(ii) = string(dateFmt);
        T_AxyFiles.timeFmt(ii) = string(timeFmt);
        T_AxyFiles.dataLines(ii) = {string(dataLines)};
        
        if hasError
            hold on;
        end
        
        [status,cmdout] = system(sprintf('wc -l "%s"',fname)); % unix
        if status == 0
            cmdParts = strsplit(strtrim(cmdout));
            T_AxyFiles.rows(ii) = str2double(cmdParts(1));
        else
            T_AxyFiles.rows(ii) = NaN;
        end
        if ~isnan(T_AxyFiles.Fs(ii))
            T_AxyFiles.days(ii) = floor(T_AxyFiles.rows(ii) / (86400*T_AxyFiles.Fs(ii)));
        else
            T_AxyFiles.days(ii) = NaN;
        end
        fprintf("%03d/%i (%1.1f%%) %s - %i days\n",ii,size(files,1),100*ii/size(files,1),T_AxyFiles.filename(ii),T_AxyFiles.days(ii));
    end
    warning ('on','all');
    do = 0;
    save("T_AxyFiles",'T_AxyFiles','rootDir');
    t = T_AxyFiles; % backup variable
    toc;
end
fprintf("%i/%i errors, %i/%i have headers\n",sum(T_AxyFiles.hasError),size(T_AxyFiles,1),sum(T_AxyFiles.hasHeader),size(T_AxyFiles,1));

%% find duplicates
sameList = [];
sameCount = 0;
clc
for ii = 1:size(T_AxyFiles)
    tryIdx = find(T_AxyFiles.rows == T_AxyFiles.rows(ii) & T_AxyFiles.id ~= ii & ~ismember(ii,sameList(:)));
    for jj = 1:numel(tryIdx)
        if all((T_AxyFiles.dataLines{ii} == T_AxyFiles.dataLines{tryIdx(jj)}))
            sameCount = sameCount + 1;
            sameList(sameCount,1:2) = [ii tryIdx(jj)]; %#ok<SAGROW>
            fprintf("%s (%s)\n%s (%s)\n\n",T_AxyFiles.filename(ii),T_AxyFiles.folder(ii),...
                T_AxyFiles.filename(tryIdx(jj)),T_AxyFiles.folder(tryIdx(jj)));
        end
    end
end
fprintf("%i/%i entries the same\n",size(sameList,1),size(T_AxyFiles,1));
%%
answer = questdlg("Remove duplicates?","Duplicates","Yes","No","No");
if strcmp(answer,"Yes")
    T_AxyFiles(sameList(:,2),:) = [];
    T_AxyFiles.id(1:size(T_AxyFiles,1)) = 1:size(T_AxyFiles,1);
    fprintf("Removed %i files, now contains %i files\n",size(sameList,1),size(T_AxyFiles,1));
else
    fprintf("No files removed\n");
end

%%
close all;
rows = 1;
cols = 2;
ff(300*cols,300);
subplot(rows,cols,1);
histogram(T_AxyFiles.days);
title('Days per session');
grid on;

subplot(rows,cols,2);
histogram(year(T_AxyFiles.startDate),min(year(T_AxyFiles.startDate))-0.5:max(year(T_AxyFiles.startDate))+0.5);
title('Sessions per year');
xtickangle(30);
grid on;

fprintf("%i cumulative days\n",sum(T_AxyFiles.days));