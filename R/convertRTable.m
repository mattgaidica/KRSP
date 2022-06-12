seasonLabels = {'Winter','Spring','Summer','Autumn'};
mastLabels = {'Non-mast','Mast'};

filename = '/Users/matt/Dropbox (Personal)/Science/KRSP Axy Manuscript/QB_pairwise_season-mast.csv';
T = readtable(filename);
VariableNames = {'Comparison','Estimate','SE','df','t','p'};
data = cell(size(T,1),numel(VariableNames));
Tnew = cell2table(data,'VariableNames',VariableNames);

for iRow = 1:size(T,1)
    Tnew{iRow,1} = {sprintf('%s %s - %s %s',seasonLabels{T{iRow,1}},mastLabels{T{iRow,2}+1},seasonLabels{T{iRow,4}},mastLabels{T{iRow,5}+1})};
    Tnew{iRow,2} = {T{iRow,6}};
    Tnew{iRow,3} = {T{iRow,7}};
    Tnew{iRow,4} = {T{iRow,8}};
    Tnew{iRow,5} = {T{iRow,9}};
    Tnew{iRow,6} = {T{iRow,10}};
end
writetable(Tnew,[filename,'.xlsx']);

%%
filename = '/Users/matt/Dropbox (Personal)/Science/KRSP Axy Manuscript/RI_pairwise_season-mast.csv';
T = readtable(filename);
VariableNames = {'Comparison','Estimate','SE','df','t','p'};
data = cell(size(T,1),numel(VariableNames));
Tnew = cell2table(data,'VariableNames',VariableNames);

for iRow = 1:size(T,1)
    Tnew{iRow,1} = {sprintf('%s %s - %s %s',seasonLabels{T{iRow,1}},mastLabels{T{iRow,2}+1},seasonLabels{T{iRow,4}},mastLabels{T{iRow,5}+1})};
    Tnew{iRow,2} = {T{iRow,6}};
    Tnew{iRow,3} = {T{iRow,7}};
    Tnew{iRow,4} = {T{iRow,8}};
    Tnew{iRow,5} = {T{iRow,9}};
    Tnew{iRow,6} = {T{iRow,10}};
end
writetable(Tnew,[filename,'.xlsx']);

%% !! skip all except mast-non-mast comparisons
filename = '/Users/matt/Dropbox (Personal)/Science/KRSP Axy Manuscript/QBnest_pairwise_season-mast.csv';
T = readtable(filename);
VariableNames = {'Comparison','Estimate','SE','df','t','p'};
data = cell(size(T,1),numel(VariableNames));
Tnew = cell2table(data,'VariableNames',VariableNames);

for iRow = 1:size(T,1)
    Tnew{iRow,1} = {sprintf('%s %s - %s %s',seasonLabels{T{iRow,1}},mastLabels{T{iRow,2}+1},seasonLabels{T{iRow,4}},mastLabels{T{iRow,5}+1})};
    Tnew{iRow,2} = {T{iRow,6}};
    Tnew{iRow,3} = {T{iRow,7}};
    Tnew{iRow,4} = {T{iRow,8}};
    Tnew{iRow,5} = {T{iRow,9}};
    Tnew{iRow,6} = {T{iRow,10}};
end
writetable(Tnew,[filename,'.xlsx']);