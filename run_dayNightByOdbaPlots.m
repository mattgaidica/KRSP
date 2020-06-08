if ismac
    filespath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
else
    filespath = 'C:\Users\mgaidica\Box Sync\KRSP Axy Data\Temp';
end
files = dir(fullfile(filespath,'*.csv.mat'));
doFig = true;
doWrite = true;
for iFile = 2:numel(files)
    loadfile = fullfile(filespath,files(iFile).name);
    disp(loadfile);
    find_dayNight(loadfile,doFig,doWrite);
end