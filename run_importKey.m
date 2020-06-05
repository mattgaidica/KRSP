% sqkey = readtable('/Users/matt/Documents/Data/KRSP/SQRaxy_key_mid_den.csv');
squirrels = unique(sqkey.Squirrel);
sqCount = [];
for iS = 1:numel(squirrels)
    sqCount(iS) = sum(strcmp(sqkey.Squirrel,squirrels{iS}));
end
close all
ff(600,500);
histogram(sqCount);
title('multiple recordings');
xlabel('# recordings');
ylabel('# squirrels');