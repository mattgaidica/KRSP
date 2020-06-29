copyFiles = false;
if copyFiles
    filespath = '/Volumes/Seagate Expansion Drive 1/2019 AXY data/Axy Data';
    files = dir2(filespath,'*.csv','-r');
    
    savePath = '/Volumes/Seagate Expansion Drive 1/2019 AXY data/2019 Exported Axy';
    
    for iFile = 1:numel(files)
        newName = strrep(files(iFile).name(1:end-3),' ','');
        newName = strrep(newName,filesep,'_');
        newName = strrep(newName,'.','');
        newName = strrep(newName,'#','');
        newName = [newName,'.csv'];
        disp(newName);
        copyfile(fullfile(filespath,files(iFile).name),fullfile(savePath,newName));
    end
end

axy19 = '/Volumes/Seagate Expansion Drive 1/2019 AXY data/2019 Exported Axy';
files19 = dir(fullfile(axy19,'*.csv'));
T=readtable('/Users/matt/Box Sync/KRSP Axy Data/Dantzer Lab Axy Data/2019 AXY data/2019 Playback Experiment-Data_16Oct_axyLogSheet_mod.csv');
for iFile = 1:numel(files19)
    strs = strsplit(files19(iFile).name,'_')
    axy = str2num(strrep(strs{1},'BJD',''));
    if ~isempty(str2num(strs{2})) && ~isempty(axy)
        matchid = find(T.SquirrelID == str2double(strs{2}) & T.Axy_ == axy);
        if numel(matchid) == 1
            T.filename{matchid} = files19(iFile).name;
        end
    end
end
% manual updates
if size(T,1) ~= 68
    error('table size changed, check manual entries!');
end
T.filename(17) = {'BJD18_22167_APRIL2019_BJD18_1.csv'}; % two dates
T.filename(48) = {'BJD10_22290_May-July2019_BJD10_1.csv'}; % two dates
T.filename(60) = {'BJD18_22167_JUNE2019_BJD18_1.csv'}; % two dates
T.filename(67) = {'BJD10_22290_October2019_October2019_1.csv'}; % two dates
T.filename(37) = {'BJD11_21922_MAY2019_BJD11_1.csv'}; % three CSVs, this is largest
T.filename(24) = {'BJD14_20033_APRILMAY2019_BJD14_2.csv'}; % three CSVs, this is largest
T.filename(45) = {'BJD14_23866_MAY-JUNE_BJD14(ATTEMPT2)_1.csv'}; % two CSVs, this is largest
T.filename(16) = {'BJD22_21944_APRIL2019_BJD22_1.csv'}; % two CSVs, this is largest
T.filename(18) = {'BJD5_22803_APRILMAY2019_BJD5_1.csv'}; % two CSVs, this is largest

% missing homes
missingHomes = {};
missCount = 0;
for iFile = 1:numel(files19)
    if ~any(strcmp(T.filename,files19(iFile).name))
        missCount = missCount + 1;
        missingHomes{missCount,1} = files19(iFile).name;
    end
end
% accounted for
missingHomes(strcmp('BJD13_20488_May2019_BJD13_1.csv',missingHomes)) = []; % mult files
missingHomes(strcmp('BJD11_21922_MAY2019_BJD11_2.csv',missingHomes)) = []; % mult files
missingHomes(strcmp('BJD11_21922_MAY2019_BJD11_3.csv',missingHomes)) = []; % mult files
missingHomes(strcmp('BJD13_23359_APRILMAY2019_BJD13_1.csv',missingHomes)) = []; % mult files
missingHomes(strcmp('BJD13_BlankData_BJD13_Blank_1.csv',missingHomes)) = []; % blank
missingHomes(strcmp('BJD14_20033_APRILMAY2019_BJD14_1.csv',missingHomes)) = []; % mult files
missingHomes(strcmp('BJD14_20033_APRILMAY2019_BJD14_3.csv',missingHomes)) = []; % mult files
missingHomes(strcmp('BJD14_23866_MAY-JUNE_BJD14_1.csv',missingHomes)) = []; % mult files
missingHomes(strcmp('BJD15_21944_June-July2019_BJD15_1.csv',missingHomes)) = []; % mult files
missingHomes(strcmp('BJD22_21944_APRIL2019_BJD22_2.csv',missingHomes)) = []; % mult files
missingHomes(strcmp('BJD5_22803_APRILMAY2019_BJD5_2.csv',missingHomes)) = []; % mult files
missingHomes(strcmp('BJD6_21465_MAY-JULY2019_BJD06_1.csv',missingHomes)) = []; % mult files

missingHomes

writetable(T,'2019modexport.csv');