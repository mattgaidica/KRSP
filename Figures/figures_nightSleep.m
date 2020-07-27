if do
    sqkey = readtable('sqkey');
    filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
    Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
    Tss_doys = day(Tss.sunrise,'dayofyear');
    sq_ids = [];
    sq_years = [];
    sq_dayLength = [];
    iRow = 0;
    corr_dist = 12;
    corr_win = 1;
    lag_tally = [];
    rs_tally = [];
    squirrelIds = [];
    sq_months = [];
    sq_sex = [];
    sq_doys = [];
    sq_days = NaT;
    invalidCount = 0;
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq}) % && ~any(ismember(sqkey.year(iSq),[2014,2019])) % ~(strcmp(sqkey.source{iSq},'ES') &&
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            T = detect_sleepWake(T);
            if isValidT(T,false)
                if numel(T.datetime) > (1440 * 3) % more than two days
                    corr_win = 60*8;
%                     awake_bin = T.odba < 1 & T.odba > 0.25; % T.awake;
                    awake_bin = T.odba > 1; % T.awake;
                    awake_sum = movsum(awake_bin,corr_win,'omitnan');
                    asleep_bin = abs(~T.awake);
                    asleep_sum = movsum(asleep_bin,corr_win,'omitnan');
                    [r,lags] = xcorr(awake_sum,asleep_sum,1440);
                    [locs,pks] = peakseek(r-mean(r),500,0);
                    if numel(locs) == 2
                        [rho1,p1] = corr(awake_sum,lagmatrix(asleep_sum,lags(locs(1))),'rows','complete');
                        [rho2,p2] = corr(awake_sum,lagmatrix(asleep_sum,lags(locs(2))),'rows','complete');
                        iRow = iRow + 1;
                        rs_tally(iRow,:) = [rho1 rho2];
                        lag_tally(iRow,:) = lags(locs);
                        sq_ids(iRow) = iSq;
                        sq_months(iRow) = month(T.datetime(1));
                        sq_sex(iRow) = strcmp(sqkey.sex{iSq},'M');
                        sq_doys(iRow) = day(T.datetime(1),'dayofyear');
                        sq_days(iRow) = T.datetime(1);
                    else
                        invalidCount = invalidCount + 1;
                    end
                end
            end
        end
    end
%     do = false;
end
fprintf('invalid: %i\n',invalidCount);

% pos=blue=sleep affects awake
nBins = 25; % make odd
bw = 360;
binEdges = linspace(720-bw,720+bw,nBins+1);
t = linspace(720-bw,720+bw,nBins);
lns = [];
lw = 2;
lag_labels = {'Awake \rightarrow Asleep','Asleep \rightarrow Awake'};
% close all
rows = 3;
cols = 4;
titleLabels = {'winter','spring','summer','fall'};
seasons = reshape(circshift(1:12,2),[3,4])';
ff(1500,1200,2);

for iSeason = 1:4
    subplot(rows,cols,prc(cols,[1,iSeason]));
    useIds = ismember(sq_months,seasons(iSeason,:));
    awake_asleep = histcounts(abs(lag_tally(useIds,1)),binEdges);
    asleep_awake = histcounts(lag_tally(useIds,2),binEdges);
    lns(1) = plot(t,awake_asleep,'linewidth',lw);
    hold on;
    lns(2) = plot(t,asleep_awake,'linewidth',lw);
    xticks(min(binEdges):60:max(binEdges));
    xticklabels(compose('%1.0f',(xticks)/60));
    xtickangle(30);
    xlim([min(t) max(t)]);
    plot([720 720],ylim,'k:');
    xlabel('hours');
    ylabel('observations');
    set(gca,'fontsize',14);
    legend(lns,lag_labels,'location','northeast');
    legend box off;
    p = anova1(abs(lag_tally(useIds,:)),[],'off');
    title(sprintf('%s\nLag of Maximum Correlation\n(p = %1.2e)',titleLabels{iSeason},p));

    subplot(rows,cols,prc(cols,[2,iSeason]));
    boxplot(rs_tally(useIds,:),'colors',lines(2));
    xticklabels({});
    set(gca,'fontsize',14);
    ylabel('corr coeff. (r)');
    xtickangle(30);
    p = anova1(rs_tally(useIds,:),[],'off');
    title(sprintf('r Distribution (p = %1.2e)',p));
    ylim([-0.5 1.25]);
    yticks([-0.5:0.5:1]);
    if p > 0.05
        H = sigstar([1,2],NaN);
        set(H(2),'fontsize',14);
    else
        H = sigstar([1,2],p);
        set(H(2),'fontsize',22);
    end
    set(H(2),'verticalalignment','bottom');
    hold on;
    plot(xlim,[0 0],'k:');

    subplot(rows,cols,prc(cols,[3,iSeason]));
    n_awake_asleep = sum(abs(rs_tally(useIds,1)) > rs_tally(useIds,2));
    n_asleep_awake = sum(abs(rs_tally(useIds,1)) < rs_tally(useIds,2));
    bar(1:2,diag([n_awake_asleep,n_asleep_awake]),'stacked');
    xticklabels(lag_labels);
    set(gca,'fontsize',14);
    ylabel('observations');
    xtickangle(30);
    [h,p] = ttest([ones(1,n_awake_asleep)*-1,ones(1,n_asleep_awake)]);
    title(sprintf('Strongest Correlation (p = %1.2e)',p));
    if p > 0.05
        H = sigstar([1,2],NaN);
        set(H(2),'fontsize',14);
    else
        H = sigstar([1,2],p);
        set(H(2),'fontsize',22);
    end
end

%%
% really similar across all seasons! I wouldn't go down that rabbit hole
% figure;plot(sq_days,lag_tally(:,2),'k.')
% shows that our data are too sparse: relationships are driven by single
% years/studies
close all
ff(1200,400);
seasons = reshape(circshift(1:12,1),[3,4])';
y = [];
group = [];
seasonVariance = [];
for iSeason = 1:4
    thisDist = abs(lag_tally(ismember(sq_months,seasons(iSeason,:)),2));
    seasonVariance(iSeason) = std(thisDist);
    y = [y thisDist'];
    group = [group repmat(iSeason,[1,numel(thisDist)])];
    subplot(1,4,iSeason);
    n = histcounts(thisDist,binEdges);
    plot(t,n);
    hold on;
    plot([720 720],ylim,'k:');
end
figure;
[p,~,stats] = anova1(y,group);
figure;
[c,m,h,nms] = multcompare(stats);

figure;
plot(seasonVariance,'-','linewidth',3);
%%
filtIds = zeros(1,size(sq_odba,1)); % blanket filter
for ii = 1:size(sq_odba,1)
    if sum(sq_odba_max(ii,300:600)) < 300
        filtIds(ii) = 1;
    end
end

corr_x = [];
corr_y = [];
dayArr = 1:366;
nDays = 7;
iRow = 0;
all_means = NaN(366,1);
all_std = NaN(366,1);
all_awakes = NaN(366,1);
for iDoy = dayArr
    shiftArr = circshift(dayArr,1-iDoy+round(nDays/2));
    useDoys = shiftArr(1:nDays);
    mastIds = ismember(sq_doys,useDoys) & filtIds; % used for awake timing
    if sum(mastIds) >= 10 % must have 10 days
        theseOdbas = sq_odba(mastIds,:);
        theseAwakes = sq_awake(mastIds,:);
        
        zGaussDay = normalize(diff(smoothdata(mean(theseAwakes),'gaussian',240,'omitnan')));
        awakeId = find(zGaussDay(600:end) > 1,1,'first')+599; % around sunrise
        % %             awakeId = 720; % sanity check for comparion
        asleepIdAM = max([1,find(zGaussDay(1:awakeId) < -1,1,'last')]); % looking back from awakeId
        asleepIdPM = min([1440,find(zGaussDay(awakeId:end) < -1,1,'first')+awakeId-1]); % forward
        
        asleepRange = asleepIdAM:awakeId;
        awakeRange = awakeId:asleepIdPM;
        
%         sum(abs(diff(theseAwakes(:,asleepRange),1,2)),2);
        iRow = iRow + 1;
        asleepOdbas = theseOdbas(:,asleepRange);
        awakeOdbas = theseOdbas(:,awakeRange);
        temp = [];
        for ii = 1:size(asleepOdbas,1)
            theseOdba = asleepOdbas(ii,:);
            theseOdba(theseOdba > 0.2) = NaN;
            temp(ii,:) = theseOdba;
        end
        asleepOdbas = temp;
        all_means(iRow) = nanmean(nanmean(asleepOdbas));
        all_std(iRow) = std(nanmean(asleepOdbas));
        all_awakes(iRow) = mean(mean(theseAwakes(:,asleepRange)));
    end
end
close all
ff(1200,900);
subplot(311);
plot(all_means,'k.')
subplot(312);
plot(all_std,'k.')
subplot(313);
plot(all_awakes,'k.')