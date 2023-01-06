function dataArr = delimText(lineText)
    for ii = 1:size(lineText,1)
        colText = strsplit(lineText(ii),"\t"); % try tab delim
        if numel(colText) == 1 % try again
            colText = strsplit(lineText(ii),","); % try comma delim
        end
        if numel(colText) == 1 % try again
            colText = strsplit(lineText(ii)," "); % try space delim
        end
        dashIdx = find(count(colText,'-') > 1); % must be a date
        colText(dashIdx) = strrep(colText(dashIdx),'-','/');
        dataArr(ii,:) = colText; %#ok<AGROW>
    end
end