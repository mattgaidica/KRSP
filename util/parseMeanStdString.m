function [m,s] = parseMeanStdString(str)
sep = '±';

m = [];
s = [];
for ii = 1:numel(str)
    sepIdx = strfind(str{ii},sep);
    if isempty(sepIdx)
        continue;
    end
    m(ii) = str2double(strrep(strtrim(str{ii}(1:sepIdx-1)),'%',''));
    s(ii) = str2double(strrep(strtrim(str{ii}(sepIdx+1:end)),'%',''));
end
m = nanmean(m);
s = nanmean(s);