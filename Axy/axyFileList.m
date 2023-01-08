if do
    doParfor = false;
    clc;
    rootDir = '/Volumes/GAIDICASSD/KRSP/KRSP Axy Data';
    csvFiles = dir2(rootDir,'-r','*.csv');
    mbSize = 1048576; % MB
    csvFiles = struct2table(csvFiles);
%     files = struct2table(files([files.bytes] >= 4096));
%     files = struct2table(files([files.bytes] >= mbSize*30));
    
    logFiles = {'2019axy_dateIssues.csv','2020 AXY logsheet.csv','AxyLog2017.csv',...
        'SQRaxy_key_mid_den.csv'};
    fnames = string;
    jj = 0;
    rmIds = [];
    for ii = 1:size(csvFiles,1)
        [~,name,ext] = fileparts(csvFiles.name(ii));
        if ~any(strcmp(name(1),{'.','~'})) && ~ismember([name,ext],logFiles)
            jj = jj + 1;
            fnames(jj,:) = string(fullfile(csvFiles.folder(ii),[name,ext]));
        else
            fprintf("skipping fname %i %s\n",ii,[name,ext]);
            rmIds = [rmIds ii]; %#ok<AGROW> 
        end
    end
    csvFiles(rmIds,:) = [];

    tic;
    results = {};
    n = numel(fnames);
    if ~doParfor
        for ii = 17%:n
            disp(ii);
            results{ii} = getAxyHeader(fnames(ii)); %#ok<SAGROW> 
        end
    else
        D = parallel.pool.DataQueue;
        h = waitbar(0, 'Please wait ...');
        nUpdateWaitbar(n, h);
        afterEach(D, @nUpdateWaitbar);
        parfor ii = 1:n
            disp(ii);
            results{ii} = getAxyHeader(fnames(ii));
            send(D,1);
        end
        close(h);
    end
    toc;

    T_AxyFiles = table;
    warning ('off','all');
    for ii = 1:numel(results)
%         T_AxyFiles.id(ii) = ii;
        [path,name,ext] = fileparts(fnames(ii));
        T_AxyFiles.folder(ii) = string(path);
        T_AxyFiles.filename(ii) = string(name+ext);
        T_AxyFiles.bytes(ii) = csvFiles.bytes(ii);
        T_AxyFiles.Fs(ii) = results{ii}.Fs;
        T_AxyFiles.colNames(ii) = {results{ii}.colNames};
        T_AxyFiles.colMap(ii) = {results{ii}.colMap};
        T_AxyFiles.headerColumns(ii) = results{ii}.headerColumns;
        T_AxyFiles.headerLabels(ii) = {results{ii}.headerLabels};
        T_AxyFiles.hasError(ii) = results{ii}.hasError;
        T_AxyFiles.hasHeader(ii) = results{ii}.hasHeader;
        T_AxyFiles.startDate(ii) = results{ii}.startDate;
        T_AxyFiles.dateFmt(ii) = string(results{ii}.dateFmt);
        T_AxyFiles.timeFmt(ii) = string(results{ii}.timeFmt);
        T_AxyFiles.rows(ii) = results{ii}.rows;
        T_AxyFiles.days(ii) = results{ii}.days;
        T_AxyFiles.md5(ii) = results{ii}.md5;
        T_AxyFiles.dataLines(ii) = {string(results{ii}.dataLines)};
    end
    warning ('on','all');
    do = 0;
    save("T_AxyFiles",'T_AxyFiles','rootDir');
    t = T_AxyFiles; % backup variable
end
clc
fprintf("%i/%i errors, %i/%i have headers\n",sum(T_AxyFiles.hasError),size(T_AxyFiles,1),sum(T_AxyFiles.hasHeader),size(T_AxyFiles,1));

%% find duplicates
% find all files with higher-level xlsx
% % axyLogIds = searchForLogs(T_AxyFiles,rootDir);

sameList = [];
sameCount = 0;
sameMenu = {};
clc
for ii = 1:size(T_AxyFiles,1)
    sameIds = find(strcmp(T_AxyFiles.md5{ii},T_AxyFiles.md5));
    if numel(sameIds) > 1 && ~any(ismember(sameIds,sameList))
        sameCount = sameCount + 1;
        sameList = [sameList;sameIds]; %#ok<AGROW> 
        sameMenu{sameCount} = sameIds; %#ok<SAGROW> 
        for jj = 1:numel(sameIds)
            fprintf("%s - %s\n",T_AxyFiles.filename(sameIds(jj)),T_AxyFiles.folder(sameIds(jj)));
        end
        fprintf("\n");
    end
end
fprintf("%i duplicates found\n",sameCount);

%!! remove duplicates? based on what?

% fprintf("%i/%i entries the same\n",size(sameList,1),size(T_AxyFiles,1));
% if sameCount > 0
%     answer = questdlg("Remove duplicates?","Duplicates","Yes","No","No");
%     if strcmp(answer,"Yes")
%         T_AxyFiles(sameList(:,2),:) = [];
%         T_AxyFiles.id(1:size(T_AxyFiles,1)) = 1:size(T_AxyFiles,1);
%         fprintf("Removed %i files, now contains %i files\n",size(sameList,1),size(T_AxyFiles,1));
%     else
%         fprintf("No files removed\n");
%     end
% end

% do again now that dups are excluded
% % [axyLogIds,no_axyLogIds,logFileList] = searchForLogs(T_AxyFiles,rootDir);
% !! could do >> T_AxyFiles(no_axyLogIds,:) == [];
% but more files might be rm in the future
% logFileList are candidates, need to work on identifying actual logs
% (maybe these all need to have "log" in filename)

%%
close all;
rows = 1;
cols = 4;
ff(1200,300);
subplot(rows,cols,1);
histogram(T_AxyFiles.days,'facecolor','k','facealpha',1,'edgecolor','w');
xlabel('rec days');
ylabel('sessions');
title('Days per session');
grid on;

subplot(rows,cols,2);
histogram(year(T_AxyFiles.startDate),min(year(T_AxyFiles.startDate))-0.5:max(year(T_AxyFiles.startDate))+0.5,...
    'facecolor','k','facealpha',1,'edgecolor','w');
xticks(min(year(T_AxyFiles.startDate)):max(year(T_AxyFiles.startDate)));
title('Sessions per year');
ylabel('sessions');
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
saveas(gcf,'KRSP_AxyDataAnalytics.jpg');