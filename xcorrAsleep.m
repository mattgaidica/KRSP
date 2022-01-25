% per squirrel or per recording?
% per season bad because it's weighted more for longer rec sessions
sq_xcorr_l = 1440*3; % 3 days to get rhythmicity index (peak 3)
unRecs = unique(sq_ids); % these are recording ids, not squirre_ids (see sq_ids_un)

sq_xcorr = [];
sq_xcorr_odba = [];
sq_xcorr_doys = [];
sq_xcorr_yrs = [];
sq_xcorr_squirrel_ids = [];
sq_xcorr_sqrow = [];
sq_xcorr_sex = [];
sq_xcorr_is_preg = [];
sq_xcorr_ndays = [];
sq_xcorr_qb = [];
sq_xcorr_dur_days = [];
sq_xcorr_trans = [];
sq_xcorr_trans_per = [];
iCount = 0;
for iSq = 1:numel(unRecs)
    useIds = find(sq_ids == unRecs(iSq));
    if numel(sq_asleep(useIds,:)) >= sq_xcorr_l
        iCount = iCount + 1;
        theseAsleep = [];
        theseODBA = [];
        for thisId = useIds
            % can use sq_odba_z here as comparison, rescale to +/-1
            theseAsleep = [theseAsleep sq_asleep(thisId,:)]; % sq_asleep
            theseODBA = [theseODBA sq_odba(thisId,:)]; % sq_asleep
        end
        [c,sq_xcorr_lags] = xcorr(normalize(theseAsleep),sq_xcorr_l,'coeff');
        sq_xcorr(iCount,:) = c;
        [c,sq_xcorr_lags] = xcorr(normalize(theseODBA),sq_xcorr_l,'coeff');
        sq_xcorr_odba(iCount,:) = c;
        sq_xcorr_doys(iCount) = sq_doys(useIds(round(numel(useIds)/2)));
        sq_xcorr_ndays(iCount) = numel(useIds);
        sq_xcorr_yrs(iCount) = sq_years(useIds(1));
        sq_xcorr_squirrel_ids(iCount) = sq_ids_un(useIds(1));
        sq_xcorr_sqrow(iCount) = sq_sqkeyrow(useIds(1));
        sq_xcorr_sex(iCount) = sq_sex(useIds(1));
        sq_xcorr_is_preg(iCount) = sq_is_preg(useIds(1));
        sq_xcorr_qb(iCount) = sum(theseAsleep==1)/numel(theseAsleep);
        sq_xcorr_dur_days(iCount) = numel(sq_asleep(useIds,:)) / 1440;
        sq_xcorr_trans(iCount) = sum(abs(diff(theseAsleep))==1);
        sq_xcorr_trans_per(iCount) = sq_xcorr_trans(iCount) / numel(useIds);
    end
end

all_RI = [];
all_RI_odba = [];
for iSq = 1:size(sq_xcorr,1)
    [ri_locs,ri_pks] = peakseek(sq_xcorr(iSq,:),720);
    ri_idx = closest(sq_xcorr_lags(ri_locs),1440*2);
    all_RI(iSq) = ri_pks(ri_idx);
    
    [ri_locs,ri_pks] = peakseek(sq_xcorr_odba(iSq,:),720);
    ri_idx = closest(sq_xcorr_lags(ri_locs),1440*2);
    all_RI_odba(iSq) = ri_pks(ri_idx);
end

%%
% generate table: all_RI, mean season, squirrel_id, year, is_mast,
% cache_size, cone_index, longevity
if ~exist('trapping','var')
    trapping = readtable('trapping.csv');
end
all_RI_Seasons = NaN(size(all_RI));
all_RI_GridConeIndex = NaN(size(all_RI));
all_RI_Longevity = NaN(size(all_RI));
all_RI_MiddenCones = NaN(size(all_RI));
all_RI_MiddenConesDiff = NaN(size(all_RI));
all_RI_Age = NaN(size(all_RI));
all_RI_Byear =  NaN(size(all_RI));
all_RI_TrapsLife = NaN(size(all_RI));
all_RI_TrapsRec = NaN(size(all_RI));
all_RI_RecWeight = NaN(size(all_RI));
all_RI_transDayRatio = NaN(size(all_RI));
all_RI_transNight = NaN(size(all_RI));
all_RI_transDay = NaN(size(all_RI));
all_RI_day_qb = NaN(size(all_RI));
all_RI_night_qb = NaN(size(all_RI));

for iSq = 1:size(sq_xcorr,1)
    for iSeason = 1:4
        if ismember(sq_xcorr_doys(iSq),useDoys{iSeason})
            break;
        end
    end
    all_RI_Seasons(iSq) = iSeason;
    gridDataId = find(strcmp(cone_counts.grid,sqkey.grid{sq_xcorr_sqrow(iSq)}) &...
        cone_counts.year == sq_xcorr_yrs(iSq));
    if ~isempty(gridDataId)
        all_RI_GridConeIndex(iSq) = cone_counts.cone_index(gridDataId);
    end
    longevityId = find(longevity.squirrel_id == sq_xcorr_squirrel_ids(iSq));
    recMidDt = datetime(sq_xcorr_yrs(iSq),1,1) + days(sq_xcorr_doys(iSq));
    if ~isempty(longevityId)
        if year(longevity.datee(longevityId)) ~= 2021
            all_RI_Longevity(iSq) = days(longevity.datee(longevityId) - longevity.dates(longevityId));
        end
% % % %         if ismember(longevity.f2(longevityId),[4,5,11,12,22])
% % % %             all_RI_Longevity(iSq) = longevity.longevity(longevityId);
% % % %         end
        all_RI_Byear(iSq) = longevity.byear(longevityId);
        all_RI_Age(iSq) = days(recMidDt - longevity.dates(longevityId));
    end
    middenId = find(midden_cones.squirrel_id == sq_xcorr_squirrel_ids(iSq) &...
        midden_cones.year == sq_xcorr_yrs(iSq));
    if ~isempty(middenId)
        all_RI_MiddenCones(iSq) = midden_cones.cache_size_total(middenId);
        all_RI_MiddenConesDiff(iSq) = midden_cones.cache_size_new(middenId) - midden_cones.cache_size_old(middenId);
    end
    trappingIds = find(trapping.squirrel_id == sq_xcorr_squirrel_ids(iSq) & trapping.wgt > 0);
    all_RI_TrapsLife(iSq) = numel(trappingIds);
    if ~isempty(trappingIds)
        daysAway = abs(days(recMidDt - trapping.date(trappingIds)));
        daysAway(daysAway == 0) = [];
        daysAway(daysAway > floor(sq_xcorr_dur_days(iSq))) = [];
        [v,k] = sort(daysAway);
        if ~isempty(k)
            recIds = find(v <= sq_xcorr_ndays(iSq)/2);
            if ~isempty(recIds)
                all_RI_TrapsRec(iSq) = numel(recIds) / sq_xcorr_ndays(iSq);
                all_RI_RecWeight(iSq) = trapping.wgt(trappingIds(k(1)));
            end
        end
    end
    all_RI_transDayRatio(iSq) = trans_dayRatio_isq(sq_xcorr_sqrow(iSq));
    all_RI_transNight(iSq) = trans_night_isq(sq_xcorr_sqrow(iSq));
    all_RI_transDay(iSq) = trans_day_isq(sq_xcorr_sqrow(iSq));
    all_RI_day_qb(iSq) = qb_day_sq(sq_xcorr_sqrow(iSq));
    all_RI_night_qb(iSq) = qb_night_sq(sq_xcorr_sqrow(iSq));
end

mastTable = zeros(size(all_RI));
mastTable(sq_xcorr_yrs == 2014 | sq_xcorr_yrs == 2019) = 1;

% !! not sure about rmoutliers, seems that extreme values might be where
% the true relationship lies
% Outliers from midden cone counts were removed (17/335),
% defined as elements more than three scaled median absolute deviations from the median
% [a,C] = rmoutliers(all_RI_MiddenCones);
% all_RI_MiddenCones(C) = NaN;
xids = ismember(all_RI_Seasons,[1,2]);
all_RI_MiddenConesDiff(xids) = NaN;
all_RI_MiddenCones(xids) = NaN;

RITable = table;
RITable.squirrel_id = sq_xcorr_squirrel_ids';
RITable.isq = sq_xcorr_sqrow';
RITable.sex = sq_xcorr_sex';
RITable.is_preg = sq_xcorr_is_preg';
RITable.doy = sq_xcorr_doys';
RITable.days = sq_xcorr_dur_days';
RITable.year = sq_xcorr_yrs';
RITable.byear = all_RI_Byear';
RITable.is_mast = mastTable';
RITable.RI = all_RI';
RITable.RI_odba = all_RI_odba';
RITable.season = all_RI_Seasons';
RITable.longevity = normalize(all_RI_Longevity');
RITable.age = normalize(all_RI_Age');
RITable.grid_cone_index = normalize(all_RI_GridConeIndex');
RITable.midden_cones = normalize(all_RI_MiddenCones');
RITable.midden_cones_diff = normalize(all_RI_MiddenConesDiff');
RITable.traps_life = normalize(all_RI_TrapsLife');
RITable.traps_rec = normalize(all_RI_TrapsRec');
RITable.trap_wgt = normalize(all_RI_RecWeight');
RITable.qb = normalize(sq_xcorr_qb');
RITable.qb_day = normalize(all_RI_day_qb');
RITable.qb_night = normalize(all_RI_night_qb');
RITable.trans = normalize(sq_xcorr_trans');
RITable.trans_per = normalize(sq_xcorr_trans_per');
RITable.trans_day_ratio = all_RI_transDayRatio';
RITable.trans_night = normalize(all_RI_transNight');
RITable.trans_day = normalize(all_RI_transDay');

all_RI_QBNest = NaN(size(all_RI));
% add in-nest qb
for ii = 1:numel(all_RI_QBNest)
    metaId = find(RITable.isq(ii) == [overlapMeta.isq{:}]);
    if ~isempty(metaId)
        all_RI_QBNest(ii) = overlapMeta.in_asleep{metaId};
    end
end
RITable.qb_nest = all_RI_QBNest';

writetable(RITable,fullfile('R','RITable.csv'));
disp("done");

%% write growthTable for R
pregLitters = find(sqkey.rec_litterSize > 0);
litterArr = [];
RIArr = [];
QBArr = [];
growthTable = table;
iArr = 0;
iGr = 0;
for iPreg = 1:numel(pregLitters)
    RIrow = find(RITable.isq == pregLitters(iPreg));
    if ~isempty(RIrow) % could be empty for short recordings (no RI)
        iArr = iArr + 1;
        
        litterId = sqkey.rec_litterId(pregLitters(iPreg));
        juvs = find(juvenile.litter_id == litterId);
        for iJuv = 1:numel(juvs)
            growthId = find(growth.squirrel_id == juvenile.squirrel_id(juvs(iJuv)));
            if ~isempty(growthId) && ~isnan(growth.growth(growthId))
                iGr = iGr + 1;
                growthTable.squirrel_id(iGr) = RITable.squirrel_id(RIrow);
                growthTable.growth(iGr) = growth.growth(growthId);
                growthTable.qb(iGr) = RITable.qb(RIrow);
                growthTable.RI(iGr) = RITable.RI(RIrow);
                growthTable.qb(iGr) = RITable.qb(RIrow);
                growthTable.is_mast(iGr) = RITable.is_mast(RIrow);
                growthTable.doy(iGr) = RITable.doy(RIrow);
                growthTable.season(iGr) = RITable.season(RIrow);
            end
        end
        
        litterArr(iArr) = sqkey.rec_litterSize(pregLitters(iPreg));
        RIArr(iArr) = RITable.RI(RIrow);
        QBArr(iArr) = RITable.qb(RIrow);
    end
end
% close all
% ff(600,600);
% subplot(221);
% scatter(litterArr,RIArr,'filled');
% xlabel('litter size');
% ylabel('RI (raw value)');
% 
% subplot(222);
% scatter(litterArr,QBArr,'filled');
% xlabel('litter size');
% ylabel('QB (Z-score)');
% 
% subplot(223);
% scatter(growthTable.growth,growthTable.RI,'filled');
% xlabel('mean growth rate');
% ylabel('RI (raw value)');
% 
% subplot(224);
% scatter(growthTable.growth,growthTable.qb,'filled');
% xlabel('mean growth rate');
% ylabel('QB (Z-score)');

writetable(growthTable,'R/GrowthTable.csv');

%% (2x) RI plots, IN PAPER
close all
doSave = 0;

ylims = [0 0.35];
months =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',366);
h = ff(400,300);
% [sort_RI,I] = sort(all_RI);
[sort_RI,I] = sort(all_RI_odba);
sort_doy = sq_xcorr_doys(I);
sort_yr = sq_xcorr_yrs(I);

mastR = [];
nmastR = [];
for iSq = 1:size(sq_xcorr,1)
    bar(iSq,sort_RI(iSq),'EdgeColor',colors(sort_doy(iSq),:),'FaceColor',colors(sort_doy(iSq),:));
    hold on;
    if ismember(sort_yr(iSq),[2014,2019])
        ln = plot(iSq,ylims(2)-0.01,'k|','linewidth',1,'markersize',10);
    else
%         plot(iSq,0.85,'k|','linewidth',1,'markersize',10);
    end
end

xlim([1 size(sq_xcorr,1)]);
set(gca,'fontsize',14);
ylabel('Rythmicity Index (RI)');
ylim(ylims);
yticks(ylims(1):0.1:ylims(2));
text(mean(xlim),ylims(2)-0.025,{'\uparrow','Mast Year'},...
    'horizontalalignment','center','verticalalignment','top','fontsize',12);
xlabel('Recording Session (sorted by RI)');
xticks([]);
grid on;
c = colorbar('location','southoutside');
colormap(colors);
c.Limits = [0,1];
c.Ticks = linspace(0,1,12);
c.TickLabels = months;
c.TickDirection = 'out';
c.FontSize = 11;
title('Rythmicity Index by Recording Session');

set(gcf,'PaperPositionMode','auto');
if doSave
    print(gcf,'-painters','-depsc',fullfile(exportPath,'rhythmicityBySession.eps')); % required for vector lines
    saveas(gcf,fullfile(exportPath,'rhythmicityBySession.jpg'),'jpg');
    close(gcf);
end

colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
h = ff(400,300);
mast_xs = [];
mast_ys = [];
mast_cmap = [];
for iSeason = 1:4
    for iMast = 1:2
        if iMast==1 %nmast
            theseIds = ismember(sq_xcorr_doys,useDoys{iSeason}) & ~ismember(sq_xcorr_yrs,[2014,2019]);
            useOffset = -0.2;
            mast_cmap = [mast_cmap;colors(iSeason,:)];
        else %mast
            theseIds = ismember(sq_xcorr_doys,useDoys{iSeason}) & ismember(sq_xcorr_yrs,[2014,2019]);
            useOffset = 0.2;
            if sum(theseIds) == 0 % winter mast is missing
                mast_xs = [mast_xs;mast_xs+useOffset*2];
                mast_ys = [mast_ys;mast_ys];
                mast_cmap = [mast_cmap;[NaN NaN NaN]];
            else
                mast_cmap = [mast_cmap;colors(iSeason,:).^2];
            end
        end
        if sum(theseIds) ~= 0 % winter mast
            mast_xs = [mast_xs;iSeason*ones(sum(theseIds),1)+useOffset];
            mast_ys = [mast_ys;all_RI_odba(theseIds)'];
        end
    end
end
bs = beeswarm(mast_xs,mast_ys,'corral_style','omit','sort_style','square','colormap',...
    mast_cmap,'overlay_style','box','MarkerFaceAlpha',1);
ylim(ylims);
yticks(ylims(1):0.1:ylims(2));
uniq_xs = unique(mast_xs);
xticks(uniq_xs);
xticklabels({'Non-mast','Mast','Non-mast','Mast','Non-mast','Mast','Non-mast','Mast'});
xtickangle(-90);
set(gca,'fontsize',14);
ylabel('Rythmicity Index (RI)');

% get pvals from R
% % % % hold on;
% % % % barY = ylims(2) - 0.06;
% % % % pvalY = ylims(2) - 0.05;
% % % % for iSeason = 1:4
% % % %     p = anova1(mast_ys(mast_xs == uniq_xs(iSeason*2-1) | mast_xs == uniq_xs(iSeason*2)),...
% % % %         mast_xs(mast_xs == uniq_xs(iSeason*2-1) | mast_xs == uniq_xs(iSeason*2)),'off');
% % % %     if p < 0.05
% % % %         plot([uniq_xs(iSeason*2-1),uniq_xs(iSeason*2)],[barY barY],'k-','linewidth',3);
% % % %         text(mean([uniq_xs(iSeason*2-1),uniq_xs(iSeason*2)]),pvalY,sprintf('*p =\n%1.2e',p),...
% % % %             'horizontalalignment','center','verticalalignment','bottom','fontsize',11);
% % % %     elseif p > 0.99 % winter mast
% % % %         plot([uniq_xs(iSeason*2-1),uniq_xs(iSeason*2)],[barY barY],'color',repmat(0.75,[1,3]),'linewidth',3);
% % % %         text(mean([uniq_xs(iSeason*2-1),uniq_xs(iSeason*2)]),pvalY,'N/A',...
% % % %             'horizontalalignment','center','verticalalignment','bottom','fontsize',11,'color',repmat(0.75,[1,3]));
% % % %     else
% % % %         plot([uniq_xs(iSeason*2-1),uniq_xs(iSeason*2)],[barY barY],'color',repmat(0.75,[1,3]),'linewidth',3);
% % % %         text(mean([uniq_xs(iSeason*2-1),uniq_xs(iSeason*2)]),pvalY,sprintf('p =\n%1.2e',p),...
% % % %             'horizontalalignment','center','verticalalignment','bottom','fontsize',11,'color',repmat(0.75,[1,3]));
% % % %     end
% % % % end

title('Rhythmicity Index Mast vs. Non-mast Years');

% need to re-open in illustrator, had to reshapre artboard, select > same >
% fill color, delete winter mast
set(gcf,'PaperPositionMode','auto');
if doSave
    print(gcf,'-painters','-depsc',fullfile(exportPath,'rhythmicityMastYears.eps')); % required for vector lines
    saveas(gcf,fullfile(exportPath,'rhythmicityMastYears.jpg'),'jpg');
    close(gcf);
end

%% by season, zoom on ultradian rhythm SUPPLEMENT FOR RI
close all
L = size(sq_xcorr,2);
vcrit = sqrt(2)*erfinv(0.95);
lconf = -vcrit/sqrt(L);
upconf = vcrit/sqrt(L);

colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
seasonLabels = {'Winter','Spring','Summer','Autumn'};
ff(700,1100);
for iZoom = 1:2
    for iSeason = 1:4
        subplot(4,2,prc(2,[iSeason,iZoom]));
        useIds = ismember(sq_xcorr_doys,useDoys{iSeason});
        xcorrMean = mean(sq_xcorr(useIds,:));
        xcorrStd = std(sq_xcorr(useIds,:));
        plot(sq_xcorr_lags,xcorrMean,'color',colors(iSeason,:),'linewidth',1.5);
        hold on;
        
        curve1 = xcorrMean + xcorrStd;
        curve2 = xcorrMean - xcorrStd;
        plot(sq_xcorr_lags,curve1,':','color',colors(iSeason,:),'linewidth',1.5);
         plot(sq_xcorr_lags,curve2,':','color',colors(iSeason,:),'linewidth',1.5);
        
        % fille isn't working on eps export
        t2 = [sq_xcorr_lags,fliplr(sq_xcorr_lags)];
        fillArea = [curve1,fliplr(curve2)];
%         fill(t2,fillArea,colors(iSeason,:),'FaceAlpha',0.15,'EdgeColor','none');
        
        % smooth to get rid of false positives
        [locs,pks] = peakseek(smoothdata(xcorrMean,'rloess',60),60);
        
        [ri_locs,ri_pks] = peakseek(xcorrMean,720,0.1);
        r1_idx = closest(sq_xcorr_lags(ri_locs),1440);
        ri_idx = closest(sq_xcorr_lags(ri_locs),1440*2);
        if iZoom == 1
            plot(sq_xcorr_lags(ri_locs(ri_idx)),ri_pks(ri_idx),'r_','markersize',20,'linewidth',1.5);
            text(sq_xcorr_lags(ri_locs(ri_idx)),ri_pks(ri_idx)+0.05,...
                sprintf('RI = %1.2f',ri_pks(ri_idx)),'fontsize',12,...
                'horizontalalignment','center','verticalalignment','bottom');
            
            plot(sq_xcorr_lags(ri_locs(r1_idx)),ri_pks(r1_idx),'r|','markersize',20,'linewidth',1.5);
            text(sq_xcorr_lags(ri_locs(r1_idx)),ri_pks(r1_idx)+0.1,...
                sprintf('RP_1 = %1.2f',sq_xcorr_lags(ri_locs(r1_idx))/60),'fontsize',12,...
                'horizontalalignment','center','verticalalignment','bottom');
            
            % 95% confidence
            yline(lconf,'r--');
            yline(upconf,'r--');
            
            xlim([min(sq_xcorr_lags) max(sq_xcorr_lags)]);
            ylim([-0.5 1]);
            title(seasonLabels{iSeason});
            xticks(-1440*3:720:1440*3);
            xticklabels(compose('%i',xticks/60));
            yticks([-0.5 0 0.5 1]);
        elseif iZoom == 2
            midVal = 720;
            midIdx = closest(sq_xcorr_lags,midVal);
            xlim([midVal-360 midVal+360]);
            ylim([xcorrMean(midIdx)-.04 xcorrMean(midIdx)+.04]);
            usePeaks = find(sq_xcorr_lags(locs) >= min(xlim) & sq_xcorr_lags(locs) <= max(xlim));
            plot(sq_xcorr_lags(locs(usePeaks(1))),xcorrMean(locs(usePeaks(1))),'r|','markersize',20,'linewidth',1.5);
            text(sq_xcorr_lags(locs(usePeaks(1))),xcorrMean(locs(usePeaks(1)))+(diff(ylim)/20),...
            compose('RP_2 = %1.2f',mod(sq_xcorr_lags(locs(usePeaks(1)))/60,24)),'fontsize',12,...
                'horizontalalignment','center','verticalalignment','bottom');
            xticks(-1440*3:180:1440*3);
            xticklabels(compose('%i',xticks/60));
            title("Ultradian Peaks");
            yticks(ylim);
            yticklabels(compose('%1.2f',yticks));
            
            leftY = ylim;
            yyaxis right;
            ylim(leftY);
            subData = xcorrMean(round(numel(xcorrMean)/2):round(numel(xcorrMean)/2)+1440-1);
            [v,k] = min(subData);
            yticks(sort([v,xcorrMean(locs(usePeaks(1)))]));
            deltaUltraPeak = xcorrMean(locs(usePeaks(1)))-v;
            set(gca,'ycolor','k');
            yticklabels([]);
            yline(v,'k:');
            yline(xcorrMean(locs(usePeaks(1))),'k:');
            vertalign = 'bottom';
            if v - min(ylim) > .01
                vertalign = 'top';
            end
            text(max(xlim)-15,v,['\DeltaRP_2 = ',sprintf('%1.4f',deltaUltraPeak)],...
                'fontsize',12,'horizontalalignment','right','verticalalignment',vertalign);
        end
        xtickangle(60)
        set(gca,'fontsize',14);
        if iSeason == 4
            xlabel('Lag (hrs)');
        end
        if iZoom == 1
            ylabel('Autocorrelation');
        end
        grid on;
    end
end

doSave = 1;
if doSave
    print(gcf,'-painters','-depsc',fullfile(exportPath,'sleepRhythmicity.eps')); % required for vector lines
    saveas(gcf,fullfile(exportPath,'sleepRhythmicity.jpg'),'jpg');
    close(gcf);
end