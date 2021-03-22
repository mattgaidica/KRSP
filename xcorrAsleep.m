sIds = round(linspace(1,366,5));
seasonDoys = circshift(1:366,57); % centers at 171, so light is equal
% seasonDoys = circshift(1:366,57-21); % centers at 192, so temp is equal
useDoys = {seasonDoys(sIds(1):sIds(2)),seasonDoys(sIds(2):sIds(3)),...
    seasonDoys(sIds(3):sIds(4)),seasonDoys(sIds(4):sIds(5))};

% two ways to do this: per squirrel or per season
% per season bad because it's weighted more for longer rec sessions
% per squirrel should still shuffle days to remove season influence
sq_xcorr_l = 1440*3; % 3 days to get rhythmicity index (peak 3)
sq_xcorr_doys = [];
sq_xcorr_yrs = [];
unSqs = unique(sq_ids);

sq_xcorr = [];
iCount = 0;
nPerm = 10;
for iSq = 1:numel(unSqs)
    useIds = find(sq_ids == unSqs(iSq));
    if numel(sq_asleep(useIds,:)) >= sq_xcorr_l
        iCount = iCount + 1;
        theseAsleep = [];
        for thisId = useIds
            % can use sq_odba_z here as comparison, rescale to +/-1
            theseAsleep = [theseAsleep 2*(sq_asleep(thisId,:)-0.5)];
        end
        [c,sq_xcorr_lags] = xcorr(theseAsleep,sq_xcorr_l,'coeff');
        sq_xcorr(iCount,:) = c;
        sq_xcorr_doys(iCount) = sq_doys(useIds(round(numel(useIds)/2)));
        sq_xcorr_yrs(iCount) = sq_years(useIds(1));
    end
end

%% all lines, colored by season
close all

% ff(1200,600);
% plot(sq_xcorr_lags,mean(sq_xcorr),'k');
% xlabel('lag (min)');
% title('all mean');

colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',366);
ff(1200,800);
for iSq = 1:size(sq_xcorr,1)
    plot(sq_xcorr_lags,sq_xcorr(iSq,:),'color',[colors(sq_doys(iSq),:),0.2],'linewidth',1);
    hold on;
end
xlim([min(sq_xcorr_lags) max(sq_xcorr_lags)]);

%% RI by squirrel
close all

all_RI = [];
for iSq = 1:size(sq_xcorr,1)
    [ri_locs,ri_pks] = peakseek(sq_xcorr(iSq,:),720);
    ri_idx = closest(sq_xcorr_lags(ri_locs),1440*2);
    all_RI(iSq) = ri_pks(ri_idx);
end

months =  {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',366);
ff(500,900);
[sort_RI,I] = sort(all_RI);
sort_doy = sq_xcorr_doys(I);
sort_yr = sq_xcorr_yrs(I);

mastR = [];
nmastR = [];
subplot(211);
for iSq = 1:size(sq_xcorr,1)
    bar(iSq,sort_RI(iSq),'EdgeColor',colors(sort_doy(iSq),:),'FaceColor',colors(sort_doy(iSq),:));
    hold on;
    if ismember(sort_yr(iSq),[2014,2019])
        ln = plot(iSq,0.975,'k|','linewidth',1,'markersize',10);
    else
%         plot(iSq,0.85,'k|','linewidth',1,'markersize',10);
    end
end

xlim([1 size(sq_xcorr,1)]);
set(gca,'fontsize',14);
ylabel('Rythmicity Index (RI)');
yticks(0:0.2:1);
text(mean(xlim),0.95,{'\uparrow','Mast Year'},...
    'horizontalalignment','center','verticalalignment','top','fontsize',12);
xlabel('Squirrel (sorted by RI)');
xticks([]);
grid on;
c = colorbar('location','southoutside');
colormap(colors);
c.Limits = [0,1];
c.Ticks = linspace(0,1,12);
c.TickLabels = months;
c.TickDirection = 'out';
c.FontSize = 11;
title('Rythmicity Index by Individual');

colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);
subplot(212);
mast_xs = [];
mast_ys = [];
mast_cmap = [];
for iSeason = 1:4
    for iMast = 1:2
        if iMast==1 %nmast
            useIds = ismember(sq_xcorr_doys,useDoys{iSeason}) & ~ismember(sq_xcorr_yrs,[2014,2019]);
            useOffset = -0.2;
            mast_cmap = [mast_cmap;colors(iSeason,:)];
        else %mast
            useIds = ismember(sq_xcorr_doys,useDoys{iSeason}) & ismember(sq_xcorr_yrs,[2014,2019]);
            useOffset = 0.2;
            if sum(useIds) == 0 % winter mast is missing, super hacky but works
                mast_xs = [mast_xs;mast_xs+useOffset*2];
                mast_ys = [mast_ys;mast_ys];
                mast_cmap = [mast_cmap;[NaN NaN NaN]];
            else
                mast_cmap = [mast_cmap;colors(iSeason,:).^2];
            end
        end
        if sum(useIds) ~= 0 % winter mast
            mast_xs = [mast_xs;iSeason*ones(sum(useIds),1)+useOffset];
            mast_ys = [mast_ys;all_RI(useIds)'];
        end
    end
end
bs = beeswarm(mast_xs,mast_ys,'corral_style','omit','sort_style','square','colormap',...
    mast_cmap,'overlay_style','box','MarkerFaceAlpha',1);
ylim([0 1]);
uniq_xs = unique(mast_xs);
xticks(uniq_xs);
xticklabels({'Non-mast','Mast','Non-mast','Mast','Non-mast','Mast','Non-mast','Mast'});
xtickangle(-90);
set(gca,'fontsize',14);
ylabel('Rythmicity Index (RI)');

hold on;
for iSeason = 1:4
    p = anova1(mast_ys(mast_xs == uniq_xs(iSeason*2-1) | mast_xs == uniq_xs(iSeason*2)),...
        mast_xs(mast_xs == uniq_xs(iSeason*2-1) | mast_xs == uniq_xs(iSeason*2)),'off');
    if p < 0.05
        plot([uniq_xs(iSeason*2-1),uniq_xs(iSeason*2)],[0.8 0.8],'k-','linewidth',3);
        text(mean([uniq_xs(iSeason*2-1),uniq_xs(iSeason*2)]),0.83,sprintf('***%1.2e',p),...
            'horizontalalignment','center','verticalalignment','bottom','fontsize',12);
    elseif p > 0.99 % winter mast
        plot([uniq_xs(iSeason*2-1),uniq_xs(iSeason*2)],[0.8 0.8],'color',repmat(0.75,[1,3]),'linewidth',3);
        text(mean([uniq_xs(iSeason*2-1),uniq_xs(iSeason*2)]),0.83,'N/A',...
            'horizontalalignment','center','verticalalignment','bottom','fontsize',12,'color',repmat(0.75,[1,3]));
    else
        plot([uniq_xs(iSeason*2-1),uniq_xs(iSeason*2)],[0.8 0.8],'color',repmat(0.75,[1,3]),'linewidth',3);
        text(mean([uniq_xs(iSeason*2-1),uniq_xs(iSeason*2)]),0.83,sprintf('%1.2e',p),...
            'horizontalalignment','center','verticalalignment','bottom','fontsize',12,'color',repmat(0.75,[1,3]));
    end
end
title('Rhythmicity Index Mast vs. Non-mast Years');

%% by season, zoom on ultradian rhythm
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
        t2 = [sq_xcorr_lags,fliplr(sq_xcorr_lags)];
        fillArea = [curve1,fliplr(curve2)];
        fill(t2,fillArea,colors(iSeason,:),'FaceAlpha',0.15,'EdgeColor','none');
        
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
        elseif iZoom == 2 || iZoom == 3
            midVal = 720;
            midIdx = closest(sq_xcorr_lags,midVal);
            xlim([midVal-360 midVal+360]);
            ylim([xcorrMean(midIdx)-.03 xcorrMean(midIdx)+.03]);
            usePeaks = sq_xcorr_lags(locs) >= min(xlim) & sq_xcorr_lags(locs) <= max(xlim);
            plot(sq_xcorr_lags(locs(usePeaks)),xcorrMean(locs(usePeaks)),'r|','markersize',20,'linewidth',1.5);
            text(sq_xcorr_lags(locs(usePeaks)),xcorrMean(locs(usePeaks))+(diff(ylim)/20),...
                compose('RP_2 = %1.2f',mod(sq_xcorr_lags(locs(usePeaks))/60,24)),'fontsize',12,...
                'horizontalalignment','center','verticalalignment','bottom');
%             if ii == 3 % 24-hr peaks
%                 midVal = 1440;
%                 midIdx = closest(sq_xcorr_lags,midVal);
%                 xlim([midVal-400 midVal+400]);
%                 xticks(-1440*3:720:1440*3);
%                 xticklabels(compose('%i',xticks/60));
%                 ylim([xcorrMean(midIdx)-xcorrStd(midIdx)-0.1 xcorrMean(midIdx)+xcorrStd(midIdx)+0.1]);
%                 title([seasonLabels{iSeason},' (zoomed)']);
%                 yticks([]);
%             else % nap peaks
                xticks(-1440*3:180:1440*3);
                xticklabels(compose('%i',xticks/60));
                title("ultradian peaks");
                yticks(ylim);
                yticklabels(compose('%1.2f',yticks));
%             end
        end
        xtickangle(60)
        set(gca,'fontsize',14);
        if iSeason == 4
            xlabel('lag (hrs)');
        end
        if iZoom == 1
            ylabel('Autocorrelation');
        end
        grid on;
    end
end