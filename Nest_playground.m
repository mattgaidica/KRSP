postMinutes = 20; % minutes
preMinutes = 20;

if do % protect
    %     do = 0;
    nestOdba = []; %NaN(1,captureMinutes);
    nestMeta = NaN(1,3);
    rowCount = 0;
    for iRec = 1:size(sq_asleep,1)
        fprintf("iRec: %i/%i\n",iRec,size(sq_asleep,1));
        if sq_sex(iRec) == 0
            curKey = sq_sqkeyrow(iRec); % keep tracking position
            continue;
        end

        nestData = sq_nest(iRec,:);
        odbaData = sq_odba(iRec,:);
        
        diffNest = diff(nestData);
        nestStateChangeIds = find(ismember(diffNest,[-1,1]));
        for ii = 1:numel(nestStateChangeIds)
            % only use 'sustained' entries/exits (based on pre/post minutes)
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
seasonColors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
unSqRow = unique(nestMeta(:,3));
stateKey = [-1,1]; % exiting, entering
sqOdba = [];
sqGradient = [];
sqSeason = [];
maxGrad = [];
sqGradDoy = [];
for iState = 1:2
    for iSq = 1:numel(unSqRow)
        useRows = find(nestMeta(:,1) == stateKey(iState) & nestMeta(:,3) == unSqRow(iSq));
        if ~isempty(useRows)
            theseOdba = nestOdba(useRows,:);
            meanOdba = nanmedian(theseOdba);
            iSeason = seasonIndex(nestMeta(useRows(1),2));
            sqSeason(iSq) = iSeason;
            sqGradDoy(iSq) = nestMeta(useRows(1),2);
        else
            meanOdba = NaN(1,size(nestOdba,2));
            sqGradDoy(iSq) = NaN;
        end
        sqOdba(iState,iSq,:) = meanOdba;
        sqGradient(iState,iSq,:) = gradient(meanOdba);
        if iState == 1
            maxGrad(iState,iSq) = max(gradient(meanOdba));
        else
            maxGrad(iState,iSq) = min(gradient(meanOdba));
        end
    end
end

% close all
% ff(1000,600);
% for iSq = 1:size(maxGrad,2)
%     plot(maxGrad(:,iSq),'color',[seasonColors(sqSeason(iSq),:),0.5],'linewidth',1.75);
%     hold on;
%     plot(1,maxGrad(1,iSq),'.','color',[seasonColors(sqSeason(iSq),:),0.5]);
%     plot(2,maxGrad(2,iSq),'.','color',[seasonColors(sqSeason(iSq),:),0.5]);
% end
% xlim([1-.1 2+.1]);
% xticks([1,2]);
% xticklabels({'Exit','Enter'});
% ylabel('ODBA Gradient');
% grid on;

close all
h1 = ff(900,900);
nestTitles = {'Just Exited Nest','Just Entered Nest'};
rows = 2;
cols = 2;
for iState = 1:2
    figure(h1);
    for iSeason = 1:4
        subplot(rows,cols,iState);
        meanOdba = squeeze(sqOdba(iState,sqSeason==iSeason,:));
        seasonOdba = nanmean(meanOdba);
        plot(seasonOdba,'color',seasonColors(iSeason,:),'linewidth',3);
        hold on;
        
        xline(preMinutes+1);
        xlim([1 preMinutes+postMinutes+1]);
        title(nestTitles{iState});
        grid on;
        xlabel('minutes');
        set(gca,'fontsize',14);
        ylabel('ODBA');
        
        subplot(rows,cols,iState+2);
        meanGradient = squeeze(sqGradient(iState,sqSeason==iSeason,:));
        seasonGradient = nanmean(meanGradient);
        plot(seasonGradient,'color',seasonColors(iSeason,:),'linewidth',3);
        hold on;
        
        xline(preMinutes+1);
        xlim([1 preMinutes+postMinutes+1]);
        title(nestTitles{iState});
        grid on;
        xlabel('minutes');
        set(gca,'fontsize',14);
        ylabel('ODBA Gradient');
    end
%     [p,anovatab,stats] = anova1(decayRates,group,'off');
%     c = multcompare(stats,'display','off');
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