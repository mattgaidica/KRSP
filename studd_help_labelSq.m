function [C,ia,ic] = studd_help_labelSq(A)
colors = "";
for ii = 1:size(A,1)
    colors(ii) = convertCharsToStrings(char([A.LCOLOUR{ii} A.RCOLOUR{ii}]));
end

[C,ia,ic] = unique(colors);