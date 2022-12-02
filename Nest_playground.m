postMinutes = 20; % minutes
preMinutes = 20;

if do % protect
    %     do = 0;
    nestOdba = []; %NaN(1,captureMinutes);
    nestMeta = NaN(1,3);
    rowCount = 0;
    for iRec = 1:size(sq_asleep,1)
        if sq_sex(iRec) == 0
            continue;
        end
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
                nestMeta(rowCount,3) = sq_sqkeyrow(iRec);
            end
        end
    end
end
%% by squirrel
unSqRow = unique(nestMeta(:,3));
stateKey = [-1,1]; % exiting, entering
sqOdba = [];
sqSeason = [];
for iState = 1:2
    for iSq = 1:numel(unSqRow)
        useRows = find(nestMeta(:,1) == stateKey(iState) & nestMeta(:,3) == unSqRow(iSq));
        theseOdba = nestOdba(useRows,:);
        meanOdba = nanmean(theseOdba);
        sqOdba(iState,iSq,:) = meanOdba;
        iSeason = seasonIndex(nestMeta(useRows(1),2));
        sqSeason(iSq) = iSeason;
    end
end

seasonColors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
close all
ff(900,600);
nestTitles = {'Just Exited Nest','Just Entered Nest'};
stateThresh = [];
seasonOdba = [];
rows = 3;
cols = 4;
figKey = [1,2,5,6;3,4,7,8];
for iState = 1:2
    subplot(rows,cols,figKey(iState,:));
    for iSeason = 1:4
        meanOdba = squeeze(sqOdba(iState,sqSeason==iSeason,:));
        seasonOdba(iSeason,:) = mean(meanOdba);
        plot(seasonOdba(iSeason,:),'color',seasonColors(iSeason,:),'linewidth',3);
        hold on;
    end
    xline(preMinutes+1);
    xlim([1 preMinutes+postMinutes+1]);
    title(nestTitles{iState});
    grid on;
    xlabel('minutes');
end

for iSeason = 1:4
    subplot(rows,cols,8+iSeason);
    y = seasonOdba(iSeason,postMinutes-2:end)';
    x = [1:numel(y)]';
    [f,gof,output] = fit(x,y,'b*x^m');
    fcoeffs = coeffvalues(f);
    confBounds = predint(f,x,0.95);
    plot(x,f(x),'color',seasonColors(iSeason,:),'linewidth',3);
    hold on;
    plot(confBounds,'-','color',seasonColors(iSeason,:));
    ylim([0.05 0.4]);
    grid on;
    title(seasonLabels{iSeason});
    legend({sprintf("m = %1.3f",fcoeffs(2))});
    xlabel('minutes');
end

%% by state/season
seasonColors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
close all
ff(1000,400);
nestTitles = {'Just Exited Nest','Just Entered Nest'};
stateKey = [-1,1]; % exiting, entering
seasonThresh = NaN(2,4);
for iState = 1:2
    subplot(1,2,iState);
    for iSeason = 1:4
        useRows = find(nestMeta(:,1) == stateKey(iState) & ismember(nestMeta(:,2),useDoys{iSeason}));
        theseOdba = nestOdba(useRows,:);
        meanOdba = nanmean(theseOdba);
        seasonThresh(iState,iSeason) = prctile(meanOdba,33);
        
        plot(meanOdba,'color',seasonColors(iSeason,:),'linewidth',3);
        hold on;
        yline(seasonThresh(iState,iSeason),'color',seasonColors(iSeason,:));
    end
    xline(preMinutes+1);
    xlim([1 preMinutes+postMinutes+1]);
    title(nestTitles{iState});
    grid on;
end

%% !! must run^ for threshold
doyColors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',366);
requireMinutes = 10;
% close all
ff(1000,400);
nestTitles = {'Just Exited Nest','Just Entered Nest'};
stateKey = [-1,1]; % exiting, entering
threshLocs = NaN(366,1);
for iState = 1:2
    subplot(1,2,iState);
    for iDoy = 1:366
        useRows = find(nestMeta(:,1) == stateKey(iState) & nestMeta(:,2) == iDoy);
        if ~isempty(useRows)
            meanOdba = smoothdata(nanmean(nestOdba(useRows,:)),'gaussian',5);
            useThresh = seasonThresh(iState,iSeason);
            if iState == 2
                threshLoc = find(meanOdba < useThresh,1,'first');
                if ~isempty(threshLoc)
                    threshLocs(iDoy) = threshLoc;
                end
            end
            plot(meanOdba,'color',[doyColors(iDoy,:) 0.2],'linewidth',1);
            hold on;
%             drawnow;
        end
    end
    xline(preMinutes+1);
    xlim([1 preMinutes+postMinutes+1]);
    ylim([0 1.4]);
    title(nestTitles{iState});
    grid on;
end
% state 2
yyaxis right;
seasonCounts = [];
for iSeason = 1:4
    seasonCounts(iSeason,:) = histcounts(threshLocs(useDoys{iSeason}),0.5:numel(meanOdba)+0.5);
    plot(seasonCounts(iSeason,:)./sum(~isnan(threshLocs(useDoys{iSeason}))),'-','color',seasonColors(iSeason,:),'linewidth',3);
end
ylabel('frac. lowest ODBA');

%%
y = [];
group = [];
jj = 1;
for iSeason = 1:4
    theseLocs = threshLocs(useDoys{iSeason});
    for ii = 1:numel(theseLocs)
        if ~isnan(theseLocs(ii))
            y(jj) = theseLocs(ii);
            group(jj) = iSeason;
            jj = jj + 1;
        end
    end
end
[p,anovatab,stats] = anova1(y,group);
multcompare(stats);