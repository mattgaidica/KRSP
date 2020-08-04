sqkey = readtable('sqkey');

for iSq = 1:size(sqkey,1)
    if ~isempty(sqkey.filename{iSq})
        load(fullfile(loadPath,sqkey.filename{iSq}));
        if ~isValidT(T,true)
            fprintf('row %i %s\n',iSq,sqkey.filename{iSq});
        end
    end
end