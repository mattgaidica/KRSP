postMinutes = 30; % minutes
preMinutes = 10;

if do % protect
%     do = 0;
    nestOdba = []; %NaN(1,captureMinutes);
    nestMeta = NaN(1,2);
    rowCount = 0;
    for iRec = 1:size(sq_asleep,1)
        fprintf("iRec: %i/%i\n",iRec,size(sq_asleep,1));
        nestData = sq_nest(iRec,:);
        odbaData = sq_odba(iRec,:);
        
        diffNest = diff(nestData);
        nestStateChangeIds = find(ismember(diffNest,[-1,1]));
        for ii = 1:numel(nestStateChangeIds)
            if (ii == numel(nestStateChangeIds) || nestStateChangeIds(ii+1) - nestStateChangeIds(ii) >= postMinutes) && nestStateChangeIds(ii) > preMinutes
                rowCount = rowCount + 1;
                nestOdba(rowCount,:) = NaN(1,preMinutes+postMinutes+1);
                nestOdba(rowCount,1:preMinutes) = odbaData(nestStateChangeIds(ii)-preMinutes:nestStateChangeIds(ii)-1);
                lastOdbaId = min([1440,nestStateChangeIds(ii)+postMinutes]);
                postData = odbaData(nestStateChangeIds(ii):lastOdbaId);
                nestOdba(rowCount,preMinutes+1:preMinutes+1+numel(postData)-1) = postData;
                nestMeta(rowCount,1) = diffNest(nestStateChangeIds(ii));
                nestMeta(rowCount,2) = sq_doys(iRec);
            end
        end
    end
end
%%
seasonColors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
requireMinutes = 10;
close all
ff(1000,400);
nestTitles = {'Just Exited Nest','Just Entered Nest'};
stateKey = [-1,1]; % exiting, entering
for iState = 1:2
    for iSeason = 1:4
        useRows = find(nestMeta(:,1) == stateKey(iState) & ismember(nestMeta(:,2),useDoys{iSeason}));
        theseOdba = nestOdba(useRows,:);
        nanCols = zeros(size(theseOdba,2),1);
        for iCol = 1:size(theseOdba,2)
            nanCols(iCol) = sum(isnan(theseOdba(:,iCol)));
        end
        meanOdba = nanmean(theseOdba);

        subplot(1,2,iState);
        plot(meanOdba,'color',seasonColors(iSeason,:),'linewidth',3);
        hold on;
        xline(preMinutes+1);
        xlim([1 preMinutes+postMinutes+1]);
        title(nestTitles{iState});
        grid on;
    end
end

%%
doyColors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',366);
requireMinutes = 10;
close all
ff(1000,400);
nestTitles = {'Just Exited Nest','Just Entered Nest'};
stateKey = [-1,1]; % exiting, entering
for iState = 1:2
    for iDoy = 1:366
        useRows = find(nestMeta(:,1) == stateKey(iState) & nestMeta(:,2) == iDoy);
        theseOdba = nestOdba(useRows,:);
        nanCols = zeros(size(theseOdba,2),1);
        for iCol = 1:size(theseOdba,2)
            nanCols(iCol) = sum(isnan(theseOdba(:,iCol)));
        end
        meanOdba = nanmean(theseOdba);

        subplot(1,2,iState);
        plot(meanOdba,'color',[doyColors(iDoy,:) 0.2],'linewidth',1);
        hold on;
        xline(preMinutes+1);
        xlim([1 preMinutes+postMinutes+1]);
        title(nestTitles{iState});
        grid on;
    end
end