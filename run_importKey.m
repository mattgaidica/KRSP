% sqkey = readtable('/Users/matt/Documents/Data/KRSP/SQRaxy_key_mid_den.csv');
% sqkey = readtable('C:\Users\mgaidica\Documents\Data\KRSP\SQRaxy_key_mid_den.csv');
% squirrels = unique(sqkey.Squirrel_ID(~isnan(sqkey.Squirrel_ID)));

axyPath = "C:\Users\mgaidica\Documents\Data\KRSP\KRSP Axy Data\Emily Studds Axy";
ff = dir(axyPath);
folds = {ff(3:end).name};
deploys = unique(sqkey.Deployment_ID);
deploy_folds = {};
vs = [];
for ii = 1:numel(deploys)
    [v,k] = min(editDistance(deploys{ii},folds));
    vs(ii) = v;
    deploy_folds{1,ii} = folds{k};
    deploy_folds{2,ii} = deploys{ii};
end
deploy_folds{1,17} = NaN;
deploy_folds{1,28} = NaN;

filenames = cell(size(sqkey,1),1);
for ii = 1:size(sqkey,1)
    inFold = deploy_folds{1,find(strcmp(sqkey.Deployment_ID{ii},{deploy_folds{2,:}}))};
    if ~isnan(inFold)
        foldPath = fullfile(axyPath,inFold);
        files = dir(fullfile(foldPath,'*.csv'));
        squirrel = sqkey.Squirrel{ii};
        for iFile = 1:numel(files)
            if strcmp(files(iFile).name(1:numel(squirrel)),squirrel)
                filenames{ii,1} = files(iFile).name;
            end
        end
    end
end
% save table



% % for iS = 1:numel(squirrels)
% %     sqCount(iS) = sum(sqkey.Squirrel_ID == squirrels(iS));
% % end
% % close all
% % ff(600,500);
% % histogram(sqCount);
% % title('multiple recordings');
% % xlabel('# recordings');
% % ylabel('# squirrels');

clc
for iS = 1:numel(squirrels)
    sqids = find(sqkey.Squirrel_ID == squirrels(iS));
    fprintf('%i:',squirrels(iS));
    for ii = 1:numel(sqids)
        if ii > 1
            fprintf('|');
        end
        fprintf('%s',sqkey.Sex{sqids(ii)})
    end
    fprintf('\n');
end