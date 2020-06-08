% sqkey = readtable('/Users/matt/Documents/Data/KRSP/SQRaxy_key_mid_den.csv');
squirrels = unique(sqkey.Squirrel_ID(~isnan(sqkey.Squirrel_ID)));
sqCount = [];

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