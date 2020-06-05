filespath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
files = dir(fullfile(filespath,'*.csv.mat'));
doFig = false;
doWrite = true;
for iFile = 1:numel(files)
    loadfile = fullfile(filespath,files(iFile).name);
    disp(loadfile);
    find_dayNight(loadfile,doFig,doWrite);
end