loadPath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
sqkey = readtable('sqkey.txt');
strFilts = {'Awake','Asleep'};
op = 0.2;
colors = lines(5);
if do
    Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
    trans_all = [];
    trans_secs = [];
    trans_months = [];
    trans_doys = [];
    trans_type = [];
    
    trans_at = [];
    trans_to = [];
    trans_on = [];
    trans_is = [];
    trans_yr = [];
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq})
            load(fullfile(loadPath,sqkey.filename{iSq}));
            if ~sqkey.isValid(iSq) % don't use temp as a test
                disp('skipping');
                continue;
            end
        else
            disp('not loading');
            continue;
        end

% % % %         if strcmp(sqkey.sex{iSq},'M') == 1
% % % %             continue;
% % % %         end

        fprintf("%i/%i - %s\n",iSq,size(sqkey,1),sqkey.filename{iSq});
        T = detect_sleepWake2(T,60);
        Tawake = make_Tawake(T); % transition table
        
% %         if ~ismember(year(Tawake.datetime(1)),[2014,2019])
% %             continue;
% %         end

        trans_at = [trans_at secDay(Tawake.datetime)'];
        trans_to = [trans_to Tawake.awake'];
        trans_on = [trans_on day(Tawake.datetime,'dayofyear')'];
        trans_is = [trans_is repmat(iSq,size(Tawake.awake'))];
        trans_yr = [trans_yr year(Tawake.datetime)'];
        
% % % %         if true % only use the following transition time
% % % %             % at
% % % %             trans_secs = [trans_secs secDay(Tawake.datetime(1:end-1))'];
% % % %             % on day
% % % %             trans_doys = [trans_doys day(Tawake.datetime(1:end-1),'dayofyear')'];
% % % %             % you transion to state
% % % %             trans_type = [trans_type Tawake.awake(2:end)'];
% % % %             % in seconds
% % % %             trans_all = [trans_all seconds(diff(Tawake.datetime))'];
% % % %         else % compile all transition times in the future
% % % %             for ii = 1:size(Tawake,1)
% % % %                 futureIdx = find(Tawake.awake(ii:end) ~= Tawake.awake(ii)) + ii - 1;
% % % %                 % limit to transitions within 24hrs, else these matrices
% % % %                 % get too big; we don't analyze further out anyways
% % % %                 futureIdx = futureIdx(seconds(Tawake.datetime(futureIdx) - Tawake.datetime(ii)) <= 86400);
% % % %                 % what if we only take out the next 2-3, assuming those
% % % %                 % might be the false positive?
% % % %                 if ~isempty(futureIdx)
% % % %                     futureIdx = futureIdx(1:min([3,numel(futureIdx)])); % modify min([X
% % % %                 end
% % % %                 
% % % %                 % at
% % % %                 trans_secs = [trans_secs repmat(secDay(Tawake.datetime(ii)),[1,numel(futureIdx)])];
% % % %                 % on day
% % % %                 trans_doys = [trans_doys repmat(day(Tawake.datetime(ii),'dayofyear'),[1,numel(futureIdx)])];
% % % %                 % you transion to state
% % % %                 trans_type = [trans_type repmat(~Tawake.awake(ii),[1,numel(futureIdx)])];
% % % %                 % in seconds
% % % %                 trans_all = [trans_all seconds(Tawake.datetime(futureIdx) - Tawake.datetime(ii))'];        
% % % %             end
% % % %         end
    end
    do = false;
%     chime;
end

%%
% extend to show 2days of repeating data
% add sun bars back in
% only plot sunrise for p(awake) and sunset for p(asleep)

Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
Tss_months = month(Tss.sunrise);
Tss_doys = day(Tss.sunrise,'dayofyear');

close all
% % % % useMonths = [1,3;4,6;7,9;10,12];
% % % % useDoys = [1,30;70,110;170,190;260,280];
sIds = round(linspace(1,366,5));
seasonDoys = circshift(1:366,60);
% useDoys = {[340:366 1:40],70:110,160:231,232:335};
useDoys = {seasonDoys(sIds(1):sIds(2)),seasonDoys(sIds(2):sIds(3)),...
    seasonDoys(sIds(3):sIds(4)),seasonDoys(sIds(4):sIds(5))};
monthNames = {'Winter','Spring','Summer','Autumn'};
typeNames = {'asleep','awake'};
h = ff;%(1400,900,2);
% % hp = ff(750,600,2);
rows = size(useDoys,2);
cols = 2;
nS = 2;

binEdges_x = linspace(0,86400,24*4+1) - (86400/2);
binEdges_y = linspace(0,16*60*60,48+1); % hours*60*60

pMats = NaN(2,size(useDoys,2),numel(binEdges_y)-1,numel(binEdges_x)-1);
for iType = 0:1
    doSunrise = true;
    trans_secsAdj = trans_secs;
    if doSunrise % center around sunrise, wrap values if necessary
        for ii = 1:numel(trans_secs)
            if iType == 0
                trans_secsAdj(ii) = trans_secs(ii) - secDay(Tss.sunset(trans_doys(ii)));
            else
                trans_secsAdj(ii) = trans_secs(ii) - secDay(Tss.sunrise(trans_doys(ii)));
            end
            if trans_secsAdj(ii) > (86400/2)
                trans_secsAdj(ii) = -86400 + trans_secsAdj(ii);
            end
            if trans_secsAdj(ii) < (-86400/2)
                trans_secsAdj(ii) = 86400 + trans_secsAdj(ii);
            end
        end
    end

    for iiDoy = 1:size(useDoys,2)
        trans_mat = zeros(numel(binEdges_x)-1,numel(binEdges_y)-1);
        trans_norm = zeros(numel(binEdges_x)-1,1);
        for iBin = 1:numel(binEdges_x)-1
            useIds = find(trans_secsAdj >= binEdges_x(iBin) & trans_secsAdj < binEdges_x(iBin+1)...
                & ismember(trans_doys,useDoys{iiDoy})...
                & trans_type == iType);
            theseVals = trans_all(useIds);
            n = histcounts(theseVals,binEdges_y);
            trans_mat(iBin,:) = n;
            trans_norm(iBin) = sum(n);
        end
        trueProb = trans_mat' ./ trans_norm';
        
        if doStats
            nSurr = 100; % should be 1000+ for production
            trans_mat_surr = zeros(numel(binEdges_x)-1,numel(binEdges_y)-1);
            trans_norm_surr = zeros(numel(binEdges_x)-1,1);
            surrProb = zeros(nSurr,numel(binEdges_y)-1,numel(binEdges_x)-1);
            for iSurr = 1:nSurr
                trans_all_surr = circshift(trans_all,randi([0 numel(trans_all)]));
                for iBin = 1:numel(binEdges_x)-1
                    useIds = find(trans_secsAdj >= binEdges_x(iBin) & trans_secsAdj < binEdges_x(iBin+1)...
                        & ismember(trans_doys,useDoys{iiDoy})...
                        & trans_type == iType);
                    theseVals = trans_all_surr(useIds);
                    n = histcounts(theseVals,binEdges_y);
                    trans_mat_surr(iBin,:) = n;
                    trans_norm_surr(iBin) = sum(n);
                end
                surrProb(iSurr,:,:) = trans_mat_surr' ./ trans_norm_surr';
            end
            
            pMat = NaN(size(trueProb));
            for ii = 1:size(trueProb,1)
                for jj = 1:size(trueProb,2)
                    pMat(ii,jj) = 1 - sum(trueProb(ii,jj) > surrProb(:,ii,jj)) / nSurr;
                end
            end
            pMats(iType+1,iiDoy,:,:) = pMat;
        else
            % load them
            disp("You are using saved stats, may not be accurate!");
            load('predict_stats');
            pMat = squeeze(pMats(iType+1,iiDoy,:,:));
        end
        
% %         figure(hp);
% %         subplot(rows,cols,prc(cols,[iiDoy,iType+1]));
% %         imagesc(pMat);
% %         set(gca,'ydir','normal');
% %         caxis([0 0.05]);
        
        figure(h);
        trueProb_filt = imgaussfilt(trueProb,nS);
        trueProb_mod = -trueProb_filt;
        trueProb_mod(pMat < 0.05) = trueProb(pMat < 0.05);
        subplot(rows,cols,prc(cols,[iiDoy,iType+1]));
        
        imagesc(repmat(binEdges_x/3600,[1,2]),1:numel(binEdges_y),repmat(trueProb_mod,[1,2]));
        
        hold on;
        set(gca,'ydir','normal');

        % these are all by hand
        xticks(-12:3:12);
        xticklabelVals = compose('%1.0f',[-12,-6,0,6,12,-6,0,6 12]);
        xticklabelVals{5} = "±12";
        xticklabels(xticklabelVals);
        
        if iiDoy == size(useDoys,2)
            if iType == 0
                xlabel('Z_t (hrs rel. to sunset)');
            else
                xlabel('Z_t (hrs rel. to sunrise)');
            end
        end
        
%         caxis([-0.149 0.15]); % make sure 0 is non-trans, manually determined
        caxis([0 0.15]);
        
%         yticks([closest(binEdges_y/3600,0),closest(binEdges_y/3600,1),closest(binEdges_y/3600,3),...
%             closest(binEdges_y/3600,6),closest(binEdges_y/3600,12),closest(binEdges_y/3600,24)]);
%         yticklabels(compose('%1.1f',[0,1,3,6,12,24]));
        ylim([1 36]);
        yticklabels(yticks/2);
        if iType == 0
            ylabel('hours');
        end
        c = colorbar;
        c.Label.String = ['p(',typeNames{iType+1},')'];
        c.Label.FontSize = 14;
        c.Limits = [0 0.15];
        title(monthNames{iiDoy});
%         colormap(magma_trans);
        colormap(magma);

%         set(gca,'ColorScale','log');
        grid on;
        
        % only for sunrise/sunset bars
        tss_range = find(ismember(Tss_doys,useDoys{iiDoy}));
        theseSunrise = secDay(Tss.sunrise(tss_range));
        theseSunset = secDay(Tss.sunset(tss_range));
        theseDaylength = Tss.day_length(tss_range) / 3600; % convert to hours for x-axis
        % do manually, harder to be fancy
        if iType == 0 % locked to sunset
            for iSun = [-6,6,18]
                plot([iSun-mean(theseDaylength)/2,iSun],[max(ylim)-0.5,max(ylim)-0.5],...
                    'linewidth',8,'color',[1 1 1 0.8]);
            end
        else % locked to sunrise
            for iSun = [-18,-6,6]
                plot([iSun,iSun+mean(theseDaylength)/2],[max(ylim)-0.5,max(ylim)-0.5],...
                    'linewidth',8,'color',[1 1 1 0.8]);
            end
        end
        
        set(gca,'fontsize',14);
        set(gca,'TitleFontSizeMultiplier',1.25);
    end
end

if doStats
    save('predict_stats','pMats');
%     doStats = false;
end