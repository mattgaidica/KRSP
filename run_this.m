% close all
dayArr = 1:366;
nDays = 30;
asleepRange = 300:600; % 240 = 4 hours
awakeRange = 900:1400; % sunrise to end?
all_rs = NaN(2,366);
all_ps = NaN(2,366);
corr_x = [];
corr_y = [];

sqs = unique(sq_ids);
jj = 0;

for iSq = 1:numel(sqs)
    useIds = find(filtIds & sq_ids == iSq);
    if numel(useIds) >= 7 && any(ismember(sq_doys(useIds),235:275))
        a = normalize(sum(sq_asleep(useIds,asleepRange),2),'range');
%         b = sum(sq_odba_max(useIds,awakeRange),2);
        b = [];
        for jj = 1:numel(useIds)
            [v,k] = sort(sq_odba_max(useIds(jj),awakeRange));
            b(jj) = mean(v(end-29:end));
        end
        
        b = normalize(b,'range')';
        corr_x = [corr_x a'];
        corr_y = [corr_y b'];
    end
end
corr_x(corr_x == 0 | corr_x == 1) = NaN;
corr_y(corr_y == 0 | corr_y == 1) = NaN;
[r,p] = corr(corr_x',corr_y','rows','complete');
close all
ff(400,400,2);
plot(corr_x',corr_y','k.');
title(sprintf('r = %1.3f, p = %1.3f',r,p));


%%
for iDoy = 1:366
    shiftArr = circshift(dayArr,1-iDoy+round(nDays/2));
    useDoys = shiftArr(1:nDays);
    mastIds = ismember(sq_doys,useDoys) & filtIds;
    if sum(mastIds) >= 10 % must have 10 days
        theseOdbas = sq_odba(mastIds,:);
        theseAwakes = sq_awake(mastIds,:);
        corr_x = sum(theseAwakes(:,asleepRange),2);
        corr_y = sum(theseOdbas(:,awakeRange),2);
        [~,I] = rmoutliers(corr_x);
        corr_x(I) = NaN;
        [~,I] = rmoutliers(corr_y);
        corr_y(I) = NaN;
        if sum(~isnan(corr_x)) > 5 && sum(~isnan(corr_y)) > 5
            [r,p] = corr(corr_x,corr_y,'rows','complete');
            all_rs(iMast,iDoy) = r;
            all_ps(iMast,iDoy) = p*(size(theseOdbas,1)^2); % bonferroni correction
        end
    end
end
