odbaArr = NaN(366,1440);
lowArr = [];
highArr = [];
for iDoy = 1:366
    useIds = ismember(sq_doys,iDoy) & filtIds;
    if sum(useIds) < 7
        continue;
    end
    theseOdbas = sq_odba(useIds,:);
    odbaArr(iDoy,:) = mean(theseOdbas);
    tArr = theseOdbas;
    tArr(tArr < 0.2 & tArr > 0.5) = NaN;
    lowArr(iDoy,:) = nanmean(tArr);
    tArr = theseOdbas;
    tArr(tArr < 0.5) = NaN;
    highArr(iDoy,:) = nanmean(tArr);
end

close all
ff(600,500);
subplot(211);
imagesc(normalize(lowArr'));
set(gca,'ydir','normal');
colormap(magma);
caxis([0 1]);

subplot(212);
imagesc(normalize(highArr'));
set(gca,'ydir','normal');
colormap(magma);
caxis([0 1]);