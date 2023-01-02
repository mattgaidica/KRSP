function ii = useCol(findStr,colNames)
    for ii = 1:numel(colNames)
        if any(strcmp(findStr,colNames{ii}))
            break;
        end
    end
end