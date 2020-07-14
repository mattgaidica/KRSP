mean_std = zeros(366,1);
mean_mean = zeros(366,1);
threshs = linspace(0.1,2.5,10);
% threshs = linspace(0.5,2.5,10);
colors = cool(numel(threshs));
ff(1000,900);
for iThresh = 1:numel(threshs)
    for iDoy = 1:366
        useIds = find(sq_doys == iDoy);
        if ~isempty(useIds)
            stdvals = [];
            meanvals = [];
            for ii = 1:numel(useIds)
                thisOdba = sq_odba(useIds(ii),:);
                stdvals(ii,:) = std(thisOdba(thisOdba > 0 & thisOdba < threshs(iThresh)));
                meanvals(ii,:) = mean(thisOdba(thisOdba > 0 & thisOdba < threshs(iThresh)));
%                 stdvals(ii,:) = std(thisOdba(thisOdba > threshs(iThresh)));
            end
            mean_mean(iDoy) = mean(meanvals);
            mean_std(iDoy) = mean(stdvals);
        end
    end
    subplot(211);
    plot(mean_std,'-','color',colors(iThresh,:));
    title('std');
    hold on;
    subplot(212);
    plot(mean_mean,'-','color',colors(iThresh,:));
    title('mean');
    hold on;
end
legend(compose('%1.2f',threshs))