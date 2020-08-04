% % % % sq_bd = readtable('/Users/matt/Documents/MATLAB/KRSP/2019modexport.csv');
KR = readtable('/Volumes/Seagate Expansion Drive/2020 XZY data/2020 AXY logsheet.csv');
dataPath = '/Volumes/Seagate Expansion Drive/2020 XZY data';
csvPath = '/Volumes/Seagate Expansion Drive/2020 XZY data/csvfiles';
files = dir2(fullfile(dataPath,'*.csv'),'-r');

if doOnce
    for iFile = 1:numel(files)
        [~,name] = fileparts(files(iFile).name);
        useId = find(strcmp(KR.DATAFILENAME,name));
        if ~isempty(useId)
            disp(name);
            copyfile(fullfile(dataPath,files(iFile).name),fullfile(csvPath,[name,'.csv']));
        end
    end
    doOnce = false;
    % manually remove double quotes using find/replace, then...
    Tcsv = readtable('/Volumes/Seagate Expansion Drive/2020 XZY data/csvfiles/BJD.14.JO.F.10..GR.LAC.csv.fix.csv',...
        'readvariablenames',0);
    writetable(Tcsv,'/Volumes/Seagate Expansion Drive/2020 XZY data/csvfiles/BJD.14.JO.F.10..GR.LAC.csv');
    categorize_odba(csvPath,'');
    compile_odba(fullfile(csvPath,'nest'),'');
end

%% do trim
dataPath = '/Volumes/Seagate Expansion Drive/2020 XZY data/csvfiles/nest';
files = dir(fullfile(dataPath,'*.mat'));
close all
clear Ttrim;
% doTrim = false;
if doTrim
    for iFile = 1:numel(files)
        load(fullfile(dataPath,files(iFile).name)); % T, Tstat, Ttrim (if there)
        if exist('Ttrim') % !! will not overwrite
            clear Ttrim;
            continue;
        end
        disp(files(iFile).name);
        h = ff(1400,500);
        plot(T.odba);
        ylabel('ODBA');
        yyaxis right;
        plot(T.temp);
        ylabel('collar temp');
        [xs,~] = ginput(2);
        close(h);
        startDate = T.datetime(round(xs(1)));
        endDate = T.datetime(round(xs(2)));
        Ttrim = [startDate,endDate];
        save(fullfile(dataPath,files(iFile).name),'T','Tstat','Ttrim');
        clear Ttrim;
    end
    doTrim = false;
end

%% apply trim dates to original CSV files
% categorize_odba so that nest class is accurate
% compile_odba using new trimmed data
% is not smart or perfect, basically needs to work once
matFiles = dir(fullfile(dataPath,'*.mat'));
savePath = '/Volumes/Seagate Expansion Drive/2020 XZY data/csvfiles/trim';
csvPath = '/Volumes/Seagate Expansion Drive/2020 XZY data/csvfiles'; % use originals
for iFile = 1:numel(matFiles)
    load(fullfile(dataPath,matFiles(iFile).name)); % Ttrim
    originalTable = readtable(fullfile(csvPath,[matFiles(iFile).name(1:end-19),'.csv']));
    startId = closest(originalTable.Var1,Ttrim(1));
    endId = closest(originalTable.Var1,Ttrim(2));
    if endId - startId > 1000
        % valid, save
        writetable(originalTable(startId:endId,:),fullfile(savePath,[matFiles(iFile).name(1:end-14),'_trim.csv']));
        disp([matFiles(iFile).name,' trimmed']);
    else
        % invalid
        disp([matFiles(iFile).name,' has invalid range']);
    end
end

% then do:
% categorize_odba('/Volumes/Seagate Expansion Drive/2020 XZY data/csvfiles/trim','_trim');
% compile_odba('/Volumes/Seagate Expansion Drive/2020 XZY data/csvfiles/trim/nest','_nest');