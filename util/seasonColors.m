function c = seasonColors(validDoys)
c = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',366);
replaceIds = ~ismember(1:366,validDoys);
c(replaceIds,:) = repmat([1 1 1],[sum(replaceIds),1]);