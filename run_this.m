%% how many records have files? This is really the starting point
iFile = 0;
for iSq = 1:size(sqkey,1)
    if isempty(sqkey.filename{iSq})
        iFile = iFile + 1;
    end
end
fprintf("%i no file\n",iFile);

%%
iValid = 0;
for iSq = 1:size(sqkey,1)
    if isempty(sqkey.filename{iSq})
        continue;
    end
%     T = loadTStruct(iSq,sqkey,Tss);
    if sqkey.isValid(iSq)
        iValid = iValid + 1;
    end
end
fprintf("%i valid\n",iValid);