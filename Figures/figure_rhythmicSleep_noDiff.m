if do
    sqkey = readtable('sqkey');
    filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
    Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
    Tss_doys = day(Tss.sunrise,'dayofyear');
    sq_odba = [];
    sq_odba_z = [];
    sq_odba_std = [];
    sq_odba_max = [];
    sq_awake = [];
    sq_ids = [];
    sq_doys = [];
    sq_years = [];
    sq_dayLength = [];
    sq_asleep = [];
    iRow = 0;
    squirrelId = 0;
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq}) % && ~any(ismember(sqkey.year(iSq),[2014,2019])) % ~(strcmp(sqkey.source{iSq},'ES') &&
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            T = detect_sleepWake2(T,60);
            if sqkey.isValid(iSq)
                dtdoys = day(T.datetime,'dayofyear');
                undoys = unique(dtdoys);
                squirrelId = squirrelId + 1;
                for iDoy = 1:numel(undoys)
% % % %                     if any(ismember(271:279,undoys(iDoy))) && strcmp(sqkey.source{iSq},'BD')
% % % %                         continue;
% % % %                     end
                    theseDoys = find(dtdoys == undoys(iDoy));
                    if numel(theseDoys) == 1440 % require full day for now
                        sunrise = Tss.sunrise(Tss_doys == undoys(iDoy));
                        afterSunrise = Tss.sunset(Tss_doys == undoys(iDoy));
                        closestId = closest(secDay(T.datetime(theseDoys)),secDay(sunrise)); % center on sunrise
                        theseDoys = theseDoys(closestId) - 720:theseDoys(closestId) + 720-1;
                        if min(theseDoys) > 1 && max(theseDoys) < numel(T.datetime)
                            iRow = iRow + 1;
                            sq_ids(iRow) = squirrelId;
                            sq_sex(iRow) = strcmp(sqkey.sex{iSq},'M'); % 0 = Female, 1 = Male
                            sq_doys(iRow) = undoys(iDoy);
                            sq_odba(iRow,:) = T.odba(theseDoys);
                            sq_odba_z(iRow,:) = T.odba_z(theseDoys);
                            sq_odba_std(iRow,:) = T.odba_std(theseDoys);
                            sq_odba_max(iRow,:) = T.odba_max(theseDoys);
                            sq_years(iRow) = year(T.datetime(theseDoys(1)));
                            sq_dayLength(iRow) = Tss.day_length(Tss_doys == undoys(iDoy));
                            sq_awake(iRow,:) = T.awake(theseDoys);
                            sq_asleep(iRow,:) = T.asleep(theseDoys);
                        end
                    end
                end
            end
        end
    end
    do = false;
end
%%

mastTitles = {'Mast','nMast'};
years_mast = [2014,2019]; % 2014,2019
years_nmast = [2015,2016,2017]; % 2015,2016,2017, *need 2018,2020
mast_years = {years_mast;years_nmast};

% create 1:366 odba statistic
pmasts = NaN(1,366);
for iDoy = 1:366
    mastIds = ismember(sq_doys,iDoy) & ismember(sq_years,years_mast) & filtIds; % sq_sex==1 & 
    nmastIds = ismember(sq_doys,iDoy) & ismember(sq_years,years_nmast) & filtIds;
    if sum(mastIds) > 1 && sum(nmastIds) > 1
        mastOdbas = sq_odba(mastIds,:);
        dist1 = sum(mastOdbas,2);
        nmastOdbas = sq_odba(nmastIds,:);
        dist2 = sum(nmastOdbas,2);
        y = [dist1;dist2];
        group = [zeros(size(dist1));ones(size(dist2))];
        pmasts(iDoy) = anova1(y,group,'off')*(numel(dist1)+numel(dist2)).^2;
    end
end
% figure;plot(find(pmasts<.05),'r*');

% close all
dayArr = 1:366;
nDays = 30;
% % % % asleepRange = 300:400; % 240 = 4 hours
% % % % awakeRange = 900:1200; % sunrise to end?
all_rs = NaN(2,366);
all_ps = NaN(2,366);
transitions = {};
odbas = {};
mastDays = zeros(2,1);
mastSquirrels = zeros(2,1);

for iMast = 1:2
    corr_x = [];
    corr_y = [];
    mastDays(iMast) = sum(ismember(sq_years,mast_years{iMast}) & filtIds);
    mastSquirrels(iMast) = numel(unique(sq_ids(ismember(sq_years,mast_years{iMast}) & filtIds)));
    for iDoy = 1:366
        shiftArr = circshift(dayArr,1-iDoy+round(nDays/2));
        useDoys = shiftArr(1:nDays);
        mastIds = ismember(sq_doys,useDoys) & ismember(sq_years,mast_years{iMast}) & filtIds;
        allIds = ismember(sq_doys,useDoys) & filtIds; % used for awake timing
        if sum(mastIds) >= 10 % must have 10 days
            theseOdbas = sq_odba(mastIds,:);
            theseAwakes = sq_awake(mastIds,:);
            
            theseAwakesAll = sq_awake(allIds,:); % so mast vs. nmast are equal vector lengths
            zGaussDay = normalize(diff(smoothdata(mean(theseAwakesAll),'gaussian',240,'omitnan')));
            awakeId = find(zGaussDay(600:end) > 1,1,'first')+599; % around sunrise
            % %             awakeId = 720; % sanity check for comparion
            asleepIdAM = max([1,find(zGaussDay(1:awakeId) < -1,1,'last')]); % looking back from awakeId
            asleepIdPM = min([1440,find(zGaussDay(awakeId:end) < -1,1,'first')+awakeId-1]); % forward
            
            asleepRange = asleepIdAM:awakeId;
            awakeRange = awakeId:asleepIdPM;
            if iMast == 1
                testSleepThresh = median(sum(abs(diff(theseAwakes(:,asleepRange),1,2)),2));
                testODBAThresh = median(sum(theseOdbas(:,awakeRange),2));
            end
            if iDoy == 265
                transitions{iMast} = abs(diff(theseAwakes(:,asleepRange),1,2));
                odbas{iMast} = theseOdbas(:,awakeRange);
                transitionsAwakeId = awakeId;
                if iMast == 1
                    h = ff(1200,800);
                    subplot(211);
                else
                    subplot(212);
                end
                plot(mean(theseOdbas));
                ylim([0 1.2]);
                yyaxis right;
                plot(zGaussDay);
                ylim([-4 4]);
                hold on;
                plot(smoothdata(normalize(mean(abs(diff(theseAwakes,1,2)))),'gaussian',50),'k-');
                plot([720,720],ylim,'y-');
                plot([awakeId,awakeId],ylim,'r-');
                plot([asleepIdAM,asleepIdAM],ylim,'g:');
                plot([asleepIdPM,asleepIdPM],ylim,'m-');
                plot(xlim,[1 1],'k:');
                plot(xlim,[-1 -1],'k:');
                xlim([-5 1446]);
                legend({'Mean ODBA','Diff ODBA','Sleep Transitions','Sunrise','"Awake"','"Asleep AM"','Asleep PM"'},...
                    'location','northwest');
                title(sprintf('%s - DOY = %i, Window = +/-%i days, n = %i days, n = %i squirrels',...
                    mastTitles{iMast},iDoy,nDays/2,size(theseAwakes,1),mastSquirrels(iMast)));
                set(gca,'fontsize',14);
                xticks(linspace(1,1440,7));
                xticklabels(linspace(-720,720,7)/60);
                xlabel(sprintf('Hours from Sunrise (%s)',datestr(Tss.sunrise(iDoy),'HH:MM AM')));
            end
            
            % this is just how many times an animal transitions
            corr_x = sum(abs(diff(theseAwakes(:,asleepRange),1,2)),2);
            % %             corr_x(corr_x < testSleepThresh) = NaN; % TEST FOR SAME SLEEP
            % %                         corr_x = randsample(corr_x,numel(corr_x)); % Sanity: SHUFFLE GETS RID OF CORRELATIONS!
            corr_y = sum(theseOdbas(:,awakeRange),2);
            % %                         corr_y(corr_y < testODBAThresh) = NaN; % TEST FOR SAME ODBA
            % do correlation
            [~,I] = rmoutliers(corr_x);
            corr_x(I) = NaN;
            [~,I] = rmoutliers(corr_y);
            corr_y(I) = NaN;
            % %             corr_x(corr_x==0) = NaN; % remove nights with no transitions?
            if sum(~isnan(corr_x)) > 5 && sum(~isnan(corr_y)) > 5
                [r,p] = corr(corr_x,corr_y,'rows','complete');
                all_rs(iMast,iDoy) = r;
                all_ps(iMast,iDoy) = p*(size(theseOdbas,1)^2); % bonferroni correction
            end
        end
    end
end


close all
ff(700,500);
lns = [];
for iMast = 1:2
    subplot(2,1,iMast);
    plot([1 366],[0 0],'-','color',[0 0 0 0.1]);
    hold on;
    plot(all_rs(iMast,:),'k.','linewidth',2);
    plot(find(all_ps(iMast,:) < 0.01),all_rs(iMast,all_ps(iMast,:) < 0.01),'b.','markersize',15);
    plot(find(all_ps(iMast,:) < 0.001),all_rs(iMast,all_ps(iMast,:) < 0.001),'r.','markersize',7);
    ylim([-1 1]);
    ylabel('corr coeff. (r)');
    
    %     yyaxis right;
    %     plot(all_ps(iMast,:),'r.');
    %     hold on;
    %     plot([1 366],[0 0],'k:');
    %     ylim([0 0.05]);
    %     ylabel('p-value');
    
    title(sprintf('(%s) Correlation between Sleep Transitions and Daytime ODBA\nDOY = %i, Window = +/-%i days, n = %i days, n = %i squirrels',...
        mastTitles{iMast},iDoy,nDays/2,mastDays(iMast),mastSquirrels(iMast)));
    xlim([1 366]);
    set(gca,'fontsize',14);
    xlabel('Day of Year');
    
    monthDoys = linspace(1,366,13);
    monthNames = {'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
    for ii = 1:12
        plot([monthDoys(ii),monthDoys(ii)],ylim,':','color',[0 0 0 0.2]);
        text(monthDoys(ii) + mean(diff(monthDoys))/2,-.83,[monthNames{ii}],...
            'fontsize',11,'horizontalalign','center','color',repmat(0.6,[1,3]));
    end
    
    lns(1) = plot(-500,-500,'b.','markersize',25);
    lns(2) = plot(-500,-500,'r.','markersize',25);
    legend(lns,{'p < 0.01','p < 0.001'});
    legend box off;
end

%% ptest for dist using surrogates
nSurr = 200;
meanBins = mean(transitions{1});
population = transitions{1}(:);
ppos = [];
pneg = [];
parr = [];
for iBin = 1:size(transitions{1},2)
    disp(iBin);
    thisMean = meanBins(iBin);
    all_dist_means = [];
    for iSurr = 1:nSurr
        all_dist_means(iSurr) = mean(randsample(population,size(transitions{1},1)));
    end
    ppos(iBin) = sum(thisMean > all_dist_means) / nSurr;
    pneg(iBin) = sum(thisMean < all_dist_means) / nSurr;
end
parr = ppos > 0.95 | pneg < 0.05;
ff(1200,300);
plot(ppos);
hold on;
plot(pneg);
plot(find(parr),0,'r*');

%% anova
y = [];
group = [];
p_bins = [];
Tmast = table;
for iBin = 1:size(transitions{1},2)
    y = [transitions{1}(:,iBin);transitions{2}(:,iBin)];
    group = [zeros(numel(transitions{1}(:,iBin)),1);ones(numel(transitions{2}(:,iBin)),1)];
    Tmast.mast(1) = sum(transitions{1}(:,iBin));
    Tmast.mast(2) = numel(transitions{1}(:,iBin)) - sum(transitions{1}(:,iBin));
    Tmast.nmast(1) = sum(transitions{2}(:,iBin));
    Tmast.nmast(2) = numel(transitions{2}(:,iBin)) - sum(transitions{2}(:,iBin));
    % these all give similar answers
    p_bins(iBin) = anova1(y,group,'off');
%     p_bins(iBin) = ranksum(transitions{1}(:,iBin),transitions{2}(:,iBin))*numel(y).^2;
%     [~,p] = fishertest(Tmast);
%     p_bins(iBin) = p;
end
parr2 = p_bins < 0.05;

% combine stats figures
t = linspace(transitionsAwakeId-720-numel(parr),transitionsAwakeId-720,numel(parr))/60;
close all
nS = 20;
lw = 2;
ms = 10;
lns = [];
ff(700,500);
subplot(211);
mastTransitions = smoothdata(sum(transitions{1}./size(transitions{1},1)),'gaussian',nS);
lns(1) = plot(t,mastTransitions,'k','linewidth',lw);
hold on;
plot(t(parr),mastTransitions(parr),'k*','markersize',ms);
nmastTransitions = smoothdata(sum(transitions{2}./size(transitions{2},1)),'gaussian',nS);
lns(2) = plot(t,nmastTransitions,'r','linewidth',lw);
hold on;
plot(t(parr2),nmastTransitions(parr2),'r*','markersize',ms);
xlim([min(t) max(t)]);
yticklabels(compose('%1.3f',yticks));
ylabel('p(transition)');
set(gca,'fontsize',14);
xlabel('Hours from Sunrise');
legend(lns,{sprintf('Mast Years (n = %i days)',size(transitions{1},1)),...
    sprintf('non-Mast Years (n = %i days)',size(transitions{2},1))},'location','southwest');
legend box off;
title('Sleep Transition Probability for DOY 265');
text(max(xlim),max(ylim),sprintf('ANOVA on means p = %1.2e',ranksum(mean(transitions{1}),mean(transitions{2}))),...
    'verticalalignment','top','horizontalalignment','right','fontsize',12);

binEdges = 1:15:numel(parr);
x = linspace(transitionsAwakeId-720-numel(parr),transitionsAwakeId-720,numel(binEdges)-1)/60;
subplot(212);
n1 = histcounts(find(parr),binEdges);
bar(x,n1,'k','facealpha',1);
hold on;
n2 = histcounts(find(parr2),binEdges);
bar(x,-n2,'r','facealpha',1);
ylim([-max([n1,n2])-1, max([n1,n2])+1]);
yticks(min(ylim):2:max(ylim));
yticklabels(abs(yticks));
ylabel('bin count');
set(gca,'fontsize',14);
xticks(-9:-2);
xlim([min(x) max(x)]);
xlabel('Hours from Sunrise');
title('Significant Deviations (p < 0.05, 15-minute bins)');
legend({'Between Mast Years','Between Mast & non-Mast Years'},'location','southwest');
legend box off;


%% shows that variation in sleep = higher odba
% potentially vigilance, anticipation, REM sleep??
% doys_mast = 95:140;
doys_mast = 235:276;
years_mast = [2014,2019];
years_nmast = [2015,2016,2017]; % need 2018,2020
filtIds = zeros(1,size(sq_odba,1)); % blanket filter
for ii = 1:size(sq_odba,1)
    if sum(sq_odba_max(ii,300:600)) < 300
        filtIds(ii) = 1;
    end
end

mastIds = ismember(sq_doys,doys_mast) & ismember(sq_years,years_mast) & filtIds; % sq_sex == 1 &
nmastIds = ismember(sq_doys,doys_mast) & ismember(sq_years,years_nmast) & filtIds;

close all
h = ff(1200,1200);
squirrelIds = {unique(sq_ids(mastIds)),unique(sq_ids(nmastIds))};
useMastIds = {mastIds,nmastIds};
mastTitles = {'Mast','nMast'};
lowMean_ind = [];
highMean_ind = [];
lowMean_coh = [];
highMean_coh = [];
lowStd_ind = [];
highStd_ind = [];
lowStd_coh = [];
highStd_coh = [];
coh_low_dist = [];
coh_high_dist = [];
all_vs = [600 600]; % do once, 600 for mean, 3000 for mean_max
all_sqvs = {};
for iMast = 1:2
    subplot(3,2,iMast);
    colors = cool(numel(squirrelIds{iMast}));
    sqvs = [];
    sqstd = [];
    sqids = [];
    for iSq = 1:numel(squirrelIds{iMast})
        thisSquirrel = sq_ids == squirrelIds{iMast}(iSq);
        theseOdbas = sq_odba(useMastIds{iMast} & thisSquirrel,:);
        theseAwakes = sq_awake(useMastIds{iMast} & thisSquirrel,:);
        % % % %         odbasStrung = reshape(cumsum(theseOdbas'),[],1);
        plot(smooth(mean(theseOdbas),100),'color',colors(iSq,:));
        hold on;
        
        if size(theseOdbas,1) > 1
            [v,k] = sort(sum(theseOdbas(:,600:end-50),2));
            sqvs = [sqvs;v]; % do once
            sqstd = [sqstd;std(theseAwakes(k,300:600),[],2)];
            sqids = [sqids;repmat(iSq,size(v))];
            
            midId = round(size(theseOdbas,1)/2);
            lowStd_ind(size(lowStd_ind,1)+1,:) = std(theseAwakes(k(1:midId),:),[],1);
            highStd_ind(size(highStd_ind,1)+1,:) = std(theseAwakes(k(midId+1:end),:),[],1);
            lowMean_ind(size(lowMean_ind,1)+1,:) = mean(theseAwakes(k(1:midId),:),1);
            highMean_ind(size(highMean_ind,1)+1,:) = mean(theseAwakes(k(midId+1:end),:),1);
            
            if size(all_vs,2) == 2 % we built it
                lowIds = v <= all_vs(iMast);
                if sum(lowIds) > 0
                    lowStd_coh(size(lowStd_coh,1)+1,:) = std(theseAwakes(lowIds,:),[],1);
                    lowMean_coh(size(lowMean_coh,1)+1,:) = mean(theseAwakes(lowIds,:),1);
                    coh_low_dist = [coh_low_dist;repmat(iSq,[sum(lowIds),1])];
                end
                highIds = v > sqvs(iMast);
                if sum(highIds) > 0
                    highStd_coh(size(highStd_coh,1)+1,:) = std(theseAwakes(highIds,:),[],1);
                    highMean_coh(size(highMean_coh,1)+1,:) = mean(theseAwakes(highIds,:),1);
                    coh_high_dist = [coh_high_dist;repmat(iSq,[sum(highIds),1])];
                end
            end
        end
    end
    % %     ff(900,400);
    % %     histogram(coh_low_dist,0.5:iSq+0.5,'edgecolor','k','displaystyle','stairs');
    % %     hold on;
    % %     histogram(coh_high_dist,0.5:iSq+0.5,'edgecolor','r','displaystyle','stairs');
    % %     title(sprintf('%s : low: %i, high: %i',mastTitles{iMast},numel(coh_low_dist),numel(coh_high_dist)));
    f = fit(sqvs,sqstd,'poly1');
    [~,I] = rmoutliers(sqvs);
    sqvs(I) = NaN;
    [~,I] = rmoutliers(sqstd);
    sqstd(I) = NaN;
    sqstd(sqstd == 0) = NaN;
    [r,p] = corr(sqvs,sqstd,'rows','complete');
    ff(500,500);
    plot(f,sqvs,sqstd);
    title(sprintf('%s : p = %1.2e, r = %1.2f',mastTitles{iMast},p,r));
    
    figure(h);
    nS = 50;
    all_sqvs{iMast} = sqvs;
    %     all_vs(iMast) = mean(sqvs);
    ylim([0 8]);
    xlim([1 size(theseOdbas,2)]);
    title(sprintf('%s (n = %i)',mastTitles{iMast},numel(sqvs)));
    
    subplot(3,2,iMast+2);
    yyaxis left; plot(smooth(mean(lowMean_ind),nS),'-k'); hold on;
    yyaxis right; plot(smooth(mean(lowStd_ind),nS),'-.k'); ylim([0 0.5]); hold on;
    yyaxis left; plot(smooth(mean(highMean_ind),nS),'-r'); hold on;
    yyaxis right; plot(smooth(mean(highStd_ind),nS),'-.r'); ylim([0 0.5]); hold on;
    title([mastTitles{iMast},' ind Awake']);
    legend('low odba','high odba','location','northwest');
    xlim([1 size(theseOdbas,2)]);
    
    subplot(3,2,iMast+4);
    yyaxis left; plot(smooth(mean(lowMean_coh),nS),'-k'); hold on;
    yyaxis right; plot(smooth(mean(lowStd_coh),nS),'-.k'); ylim([0 0.5]); hold on;
    yyaxis left; plot(smooth(mean(highMean_coh),nS),'-r'); hold on;
    yyaxis right; plot(smooth(mean(highStd_coh),nS),'-.r'); ylim([0 0.5]); hold on;
    title([mastTitles{iMast},' coh Awake']);
    legend('low odba','high odba','location','northwest');
    xlim([1 size(theseOdbas,2)]);
end

%% plot mean odba mast/nmast
awake_mast = sq_awake(mastIds,:);
awake_nmast = sq_awake(nmastIds,:);
odba_mast = sq_odba_max(mastIds,:);
odba_nmast = sq_odba_max(nmastIds,:);

todba = 1:size(odba_mast,2);
% close all
ff(1200,1200);
subplot(321);
plot(todba,mean(odba_mast));
hold on;
plot(todba,mean(odba_nmast));
legend({'mast','nmast'});
title('mean ODBA');

tawake = 1:size(awake_mast,2);
subplot(322);
plot_distribution(tawake,awake_mast);
hold on;
plot_distribution(tawake,awake_nmast);
title('mean awake');

subplot(323);
imagesc(odba_mast);
colormap(magma);caxis([0 2.5]);
title('mast');
subplot(324);
imagesc(odba_nmast);
colormap(magma);caxis([0 2.5]);
title('nmast');

subplot(325);
imagesc(awake_mast);
colormap(magma);caxis([0 1]);
title('mast');
subplot(326);
imagesc(awake_nmast);
colormap(magma);caxis([0 1]);
title('nmast');

%% all doys
mean_odba = [];
max_odba = [];
daylength_odba = [];
mean_awake = zeros(366,1440);
for iDoy = 1:366
    useIds = find(sq_doys == iDoy);
    mean_odba(iDoy,:) = mean(sq_odba(useIds,:),1);
    if ~isempty(useIds)
        mean_awake(iDoy,:) = mean(sq_awake(useIds,:),1);
        maxvals = [];
        normvals = [];
        
        A = sq_odba((sq_doys == iDoy),:);
        allA = A(:);
        todayStd = std(allA,0,1,'omitnan');
        todayMean = mean(allA);
        for ii = 1:numel(useIds)
            thisOdba = sq_odba(useIds(ii),:);
            v = sort(thisOdba);
            maxvals(ii) = std(thisOdba);
            %             thisOdba(thisOdba > 0.5) = 0.5;
            normvals(ii,:) = thisOdba;
        end
        max_odba(iDoy) = mean(maxvals');
    end
    %     daylength_odba(iDoy) = mean_odba(iDoy) ./ mean(sq_dayLength(useIds));
end

% make shifted sunset
sunsetOffset = [];
for iDoy = 1:366
    afterSunrise = Tss.day_length(iDoy)/60/60; % minutes
    
    if afterSunrise > 12
        sunsetOffset(iDoy) = afterSunrise - 24;
    else
        sunsetOffset(iDoy) = afterSunrise;
    end
end
% figure;plot(sunsetOffset,'k--','linewidth',lw);

close all
lw = 1.5;
tWindow = 1440/2;
t = linspace(-tWindow,tWindow,size(sq_awake,2))/60;
ff(700,500);
subplot(211);
imagesc(1:366,t,mean_odba');
hold on;
plot(sunsetOffset,'w--','linewidth',lw);
plot(xlim,[0,0],'-','color',[1 1 1 0.4]);
set(gca,'ydir','normal');
xlabel('Day of Year');
ylabel('Hours from Sunrise');
c = colorbar;
c.Label.String = 'Mean ODBA';
colormap(magma);
caxis([0 1]);
title('Mean ODBA');
set(gca,'fontsize',14);

subplot(212);
imagesc(1:366,t,mean_awake');
hold on;
plot(sunsetOffset,'w--','linewidth',lw);
plot(xlim,[0,0],'-','color',[1 1 1 0.4]);
set(gca,'ydir','normal');
xlabel('Day of Year');
ylabel('Hours from Sunrise');
c = colorbar;
c.Label.String = 'p(Awake)';
colormap(magma);
title('Awake Probability');
set(gca,'fontsize',14);