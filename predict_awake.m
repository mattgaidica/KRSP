loadPath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
sqkey = readtable('sqkey.txt');
strFilts = {'Awake','Asleep'};
op = 0.2;
colors = lines(5);
if do
    Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
    Tss_doys = day(Tss.sunrise,'dayofyear');
    trans_all = [];
    trans_secs = [];
    trans_months = [];
    trans_doys = [];
    trans_type = [];
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq})
            load(fullfile(loadPath,sqkey.filename{iSq}));
            if ~isValidT(T,false) % don't use temp as a test
                disp('skipping');
                continue;
            end
        end
        disp(sqkey.filename{iSq});
        T = detect_sleepWake2(T,60);
        Tawake = make_Tawake(T);
        
        % at
        trans_secs = [trans_secs secDay(Tawake.datetime(1:end-1))'];
        % on day
        trans_doys = [trans_doys day(Tawake.datetime(1:end-1),'dayofyear')'];
        % you transion to state
        trans_type = [trans_type Tawake.awake(2:end)'];
        % in seconds
        trans_all = [trans_all seconds(diff(Tawake.datetime))'];
        
        % % % %         for iT = 2:size(Tawake,1) % last entry has no transition to reference
        % % % %             transTime = seconds(Tawake.datetime(iT) - Tawake.datetime(iT-1));
        % % % %             trans_all = [trans_all transTime];
        % % % %             trans_secs = [trans_secs secDay(Tawake.datetime(iT))];
        % % % %             trans_months = [trans_months month(Tawake.datetime(iT))];
        % % % %             trans_doys = [trans_doys day(Tawake.datetime(iT),'dayofyear')];
        % % % %             trans_type = [trans_type Tawake.awake(iT)];
        % % % %         end
    end
    do = false;
    chime;
end

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
monthNames = {'Winter','Spring','Summer','Fall'};
typeNames = {'asleep','awake'};
h = ff(750,600,2);
% % hp = ff(750,600,2);
rows = size(useDoys,2);
cols = 2;
nS = 1;
binEdges_x = linspace(0,86400,24*4+1); % was 24*4+1
binEdges_y = logspace(2.9,5,48); % was 48
pMats = NaN(2,size(useDoys,2),numel(binEdges_y)-1,numel(binEdges_x)-1);
for iType = 0:1
    for iiDoy = 1:size(useDoys,2)
        trans_mat = zeros(numel(binEdges_x)-1,numel(binEdges_y)-1);
        trans_norm = zeros(numel(binEdges_x)-1,1);
        for iBin = 1:numel(binEdges_x)-1
            useIds = find(trans_secs >= binEdges_x(iBin) & trans_secs < binEdges_x(iBin+1)...
                & ismember(trans_doys,useDoys{iiDoy})...
                & trans_type == iType);
            theseVals = trans_all(useIds);
            n = histcounts(theseVals,binEdges_y);
            trans_mat(iBin,:) = n;
            trans_norm(iBin) = sum(n);
        end
        trueProb = trans_mat' ./ trans_norm';
        
        if doStats
            nSurr = 1000; % should be 1000 for production
            trans_mat_surr = zeros(numel(binEdges_x)-1,numel(binEdges_y)-1);
            trans_norm_surr = zeros(numel(binEdges_x)-1,1);
            surrProb = zeros(nSurr,numel(binEdges_y)-1,numel(binEdges_x)-1);
            for iSurr = 1:nSurr
                trans_all_surr = circshift(trans_all,randi([0 numel(trans_all)]));
                for iBin = 1:numel(binEdges_x)-1
                    useIds = find(trans_secs >= binEdges_x(iBin) & trans_secs < binEdges_x(iBin+1)...
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
        imagesc(trueProb_mod);
        hold on;
        set(gca,'ydir','normal');
        xticks(linspace(min(xlim),max(xlim),8));
        xticklabels(compose('%1.0f',linspace(0,23,numel(xticks))));
        if iiDoy == size(useDoys,2)
            xlabel('hour of day');
        end
        caxis([-0.149 0.15]); % make sure 0 is non-trans
        yticklabels(compose('%1.1f',binEdges_y(yticks)/60/60));
        if iType == 0
            ylabel('hours');
        end
        c = colorbar;
        c.Label.String = ['p(',typeNames{iType+1},')'];
        c.Label.FontSize = 14;
        c.Limits = [0 0.15];
        title(monthNames{iiDoy});
        colormap(magma_trans);
        
        % only for sunrise/sunset bars
        tss_range = find(ismember(Tss_doys,useDoys{iiDoy}));
        theseSunrise = secDay(Tss.sunrise(tss_range));
        theseSunset = secDay(Tss.sunset(tss_range));
        xsc = linspace(0,86400,max(xlim)-min(xlim));
        sunrise_x = closest(xsc,mean(theseSunrise));
        sunrise_std_x = closest(xsc,mean(theseSunrise)-std(theseSunrise));
        sunset_x = closest(xsc,mean(theseSunset));
        sunset_std_x = closest(xsc,mean(theseSunset)+std(theseSunset));
        %         plot([sunrise_std_x,sunset_std_x],[max(ylim)-0.5,max(ylim)-0.5],'linewidth',8,'color',[1 1 1 0.5]); % colors(3,:)
        plot([sunrise_x,sunset_x],[max(ylim)-0.5,max(ylim)-0.5],'linewidth',8,'color',[1 1 1 0.8]);
        
        set(gca,'fontsize',14);
        set(gca,'TitleFontSizeMultiplier',1.25);
    end
end

if doStats
    save('predict_stats','pMats');
%     doStats = false;
end