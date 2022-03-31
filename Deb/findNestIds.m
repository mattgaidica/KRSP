Nests = readtable('AllNests.csv');
usefmt = '%d%q%q%q%q%q%q%q%q%q%q%q%q%q%q';
% usefmt = '%d%D%q%q%q%D%q%D%f%f%q%q%q%q%D';
Squirrel = readtable('krsp_squirrel.csv','Format',usefmt);
%%
clc
for iNest = 1:size(Nests,1)
    if isnan(Nests.squirrel_id(iNest))
        tagId = find(strcmp(Nests.Taglft{iNest},Squirrel.taglft) & strcmp(Nests.Tagrt{iNest},Squirrel.tagrt));
        if ~isempty(tagId)
            Nests.squirrel_id(iNest) = Squirrel.id(tagId);
        else
            tagId_lft = find(strcmp(Nests.Taglft{iNest},Squirrel.taglft) & strcmp('F',Squirrel.sex));
            tagId_rt = find(strcmp(Nests.Tagrt{iNest},Squirrel.tagrt) & strcmp('F',Squirrel.sex));
            fprintf('Maybe L%i or R%i\n',numel(tagId_lft),numel(tagId_rt));
        end
    end
end
writetable(Nests,'AllNests_fix_20220330.csv');