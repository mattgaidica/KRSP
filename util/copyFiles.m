function copyFiles(fromPath,toPath,pattern)

files = dir2(fromPath,pattern,'-r');
for iFile = 1:numel(files)
    disp(files(iFile).name);
    [~,name,ext] = fileparts(files(iFile).name);
    copyfile(fullfile(fromPath,files(iFile).name),...
        fullfile(toPath,[name,ext]));
end
disp(['Done copying ',num2str(numel(files)),' files.']);