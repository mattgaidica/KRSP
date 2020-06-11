function files = dirx(files,excludeString)

keepFiles = [];
for iFile = 1:numel(files)
    if isempty(strfind(files(iFile).name,excludeString))
        keepFiles = [keepFiles;iFile];
    end
end
files = files(keepFiles);