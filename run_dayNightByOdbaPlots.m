if ismac
    filespath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
else
    filespath = 'C:\Users\mgaidica\Box Sync\KRSP Axy Data\Temp';
end
files = dirx(dir(fullfile(filespath,'*.mat')),'meta');
doFig = false;
doWrite = true;
for iFile = 1:numel(files)
    loadfile = fullfile(filespath,files(iFile).name);
    disp(loadfile(end-30:end));
    find_dayNight(loadfile,doFig,doWrite);
end