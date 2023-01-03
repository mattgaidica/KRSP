if do
    tic;
    clc;
    rootDir = '/Volumes/GAIDICASSD/KRSP/KRSP Axy Data';
    files = dir2(rootDir,'-r','*.csv');
    mbSize = 1048576; % MB
    files = struct2table(files([files.bytes] >= mbSize*30));

    fnames = string;
    for ii = 1:size(files,1)
        [~,name,ext] = fileparts(files.name(ii));
        fnames(ii,:) = string(fullfile(files.folder(ii),[name,ext]));
    end

    tic;
    results = {};
    n = numel(fnames);
    D = parallel.pool.DataQueue;
    h = waitbar(0, 'Please wait ...');
    nUpdateWaitbar(n, h);
    afterEach(D, @nUpdateWaitbar);
    parfor ii = 1:n
        results{ii} = getAxyHeader(fnames(ii));
        send(D,1);
    end
    close(h);
    toc

    T_AxyFiles = table;
    warning ('off','all');
    for ii = 1:numel(results)
        T_AxyFiles.id(ii) = ii;
        [path,name,ext] = fileparts(fnames(ii));
        T_AxyFiles.folder(ii) = string(path);
        T_AxyFiles.filename(ii) = string(name+ext);
        T_AxyFiles.mb(ii) = round(files.bytes(ii) / mbSize);
        T_AxyFiles.Fs(ii) = results{ii}.Fs;
        T_AxyFiles.colNames(ii) = {results{ii}.colNames};
        T_AxyFiles.colMap(ii) = {results{ii}.colMap};
        T_AxyFiles.hasError(ii) = results{ii}.hasError;
        T_AxyFiles.hasHeader(ii) = results{ii}.hasHeader;
        T_AxyFiles.startDate(ii) = results{ii}.startDate;
        T_AxyFiles.dateFmt(ii) = string(results{ii}.dateFmt);
        T_AxyFiles.timeFmt(ii) = string(results{ii}.timeFmt);
        T_AxyFiles.rows(ii) = results{ii}.rows;
        T_AxyFiles.days(ii) = results{ii}.days;
        T_AxyFiles.dataLines(ii) = {string(results{ii}.dataLines)};
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
% find all files with higher-level xlsx
axyLogIds = searchForLogs(T_AxyFiles,rootDir);

sameList = [];
sameCount = 0;
clc
for ii = 1:size(T_AxyFiles)
    tryIdx = find(T_AxyFiles.rows == T_AxyFiles.rows(ii) & T_AxyFiles.id ~= ii & ~ismember(ii,sameList(:)));
    for jj = 1:numel(tryIdx)
        if all((T_AxyFiles.dataLines{ii} == T_AxyFiles.dataLines{tryIdx(jj)}))
            sameCount = sameCount + 1;
            % second element will be excluded below, keep files w logs
            if ismember(ii,axyLogIds)
                sameList(sameCount,1:2) = [ii tryIdx(jj)]; %#ok<SAGROW>
            else
                sameList(sameCount,1:2) = [tryIdx(jj) ii]; %#ok<SAGROW>
            end
            fprintf("%s (%s)\n%s (%s)\n\n",T_AxyFiles.filename(ii),T_AxyFiles.folder(ii),...
                T_AxyFiles.filename(tryIdx(jj)),T_AxyFiles.folder(tryIdx(jj)));
        end
    end
end
fprintf("%i/%i entries the same\n",size(sameList,1),size(T_AxyFiles,1));
if sameCount > 0
    answer = questdlg("Remove duplicates?","Duplicates","Yes","No","No");
    if strcmp(answer,"Yes")
        T_AxyFiles(sameList(:,2),:) = [];
        T_AxyFiles.id(1:size(T_AxyFiles,1)) = 1:size(T_AxyFiles,1);
        fprintf("Removed %i files, now contains %i files\n",size(sameList,1),size(T_AxyFiles,1));
    else
        fprintf("No files removed\n");
    end
end

% do again now that dups are excluded
[axyLogIds,no_axyLogIds] = searchForLogs(T_AxyFiles,rootDir);
% !! could do >> T_AxyFiles(no_axyLogIds,:) == [];
% but more files might be rm in the future

%%
close all;
rows = 1;
cols = 4;
ff(1200,300);
subplot(rows,cols,1);
histogram(T_AxyFiles.days,'facecolor','k','facealpha',1,'edgecolor','w');
xlabel('rec days');
title('Days per session');
grid on;

subplot(rows,cols,2);
histogram(year(T_AxyFiles.startDate),min(year(T_AxyFiles.startDate))-0.5:max(year(T_AxyFiles.startDate))+0.5,...
    'facecolor','k','facealpha',1,'edgecolor','w');
xticks(min(year(T_AxyFiles.startDate)):max(year(T_AxyFiles.startDate)));
title('Sessions per year');
xtickangle(30);
grid on;

subplot(rows,cols,3:4)
recDays = zeros(size(T_AxyFiles,1),366);
[v,k] = sort(day(T_AxyFiles.startDate,'dayofyear'));
for ii = 1:size(T_AxyFiles,1)
    theseDays = zeros(1,366);
    theseDays(1:T_AxyFiles.days(k(ii))) = 1;
    theseDays = circshift(theseDays,v(ii));
    recDays(ii,:) = theseDays;
end
imagesc(recDays);
set(gca,'ydir','normal');
colormap("gray");
title('Data Coverage')
xlabel('DOY');
xticks([1,366]);
ylabel('Rec session');

fprintf("%i cumulative days\n",sum(T_AxyFiles.days));