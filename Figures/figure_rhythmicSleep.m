if do
    sqkey = readtable('sqkey');
    filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
    Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
    Tss_doys = day(Tss.sunrise,'dayofyear');
    awake_sunrise = [];
    awake_sunset = [];
    sq_odba = [];
    sq_odba_max = [];
    sq_ids = [];
    sq_doys = [];
    squirrelId = 0;
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq})
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            T = detect_sleepWake(T,2);
            awakeIds = find(diff(T.awake)==1)+1;
            if numel(awakeIds) > 2
                squirrelId = squirrelId + 1;
                awakeDts = T.datetime(awakeIds);
                dtdoys = day(T.datetime,'dayofyear');
                undoys = unique(dtdoys);
                
                for iDoy = 1:numel(undoys)
                    % get sunrise/sunset
                    theseDoys = find(dtdoys == undoys(iDoy));
                    if numel(theseDoys) == 1440 % require full day for now
                        % roughly, using 2016 right now
                        sunrise = Tss.sunrise(Tss_doys == undoys(iDoy));
                        sunset = Tss.sunset(Tss_doys == undoys(iDoy));
                        % do I want ODBA for NEXT DAY for sunset values?
                        awake_sunset = [awake_sunset;seconds(awakeDts-sunset)];
                        awake_sunrise = [awake_sunrise;seconds(awakeDts-sunrise)];
                        
                        sq_odba = [sq_odba;repmat(sum(T.odba(theseDoys)),size(awakeDts))]; % no need to normalize now
                        % %                         v = sort(T.odba(theseDoys));
                        % %                         sq_odba_max = sum(v(end-29:end)); % n to end?
                        sq_ids = [sq_ids;repmat(squirrelId,size(awakeDts))];
                        sq_doys = [sq_doys;repmat(undoys(iDoy),size(awakeDts))];
                    end
                end
            end
        end
    end
    do = false;
end

%% setup
nHours = 6;
nBins = (nHours*60*2)+1;
binEdges = linspace(-nHours*3600,nHours*3600,nBins);
nSq = numel(unique(sq_ids));
seasonDoys = round(linspace(1,366,5));
seasonLabel = {'Winter','Spring','Summer','Fall'};
close all
ff(1200,1000);
for iSeason = 1:4
    odbaHist = [];
    odbaMaxHist = [];
    sunsetHist = [];
    sunriseHist = [];
    sqDoyCount = 0;
    for iSq = 1:nSq
        unDoys = unique(sq_doys(sq_ids == iSq));
        for iDoy = 1:numel(unDoys)
            if unDoys(iDoy) >= seasonDoys(iSeason) && unDoys(iDoy) <= seasonDoys(iSeason+1)
                sqDoyCount = sqDoyCount + 1;
                odbaHist(sqDoyCount) = min(sq_odba(sq_ids == iSq & sq_doys == unDoys(iDoy)));
                % %         odbaMaxHist(sqDoyCount) = min(sq_odba_max(sq_ids == iSq & sq_doys == unDoys(iDoy)));
                counts = histcounts(awake_sunset(sq_ids == iSq & sq_doys == unDoys(iDoy)),binEdges);
                sunsetHist(sqDoyCount,:) = counts;
                counts = histcounts(awake_sunrise(sq_ids == iSq & sq_doys == unDoys(iDoy)),binEdges);
                sunriseHist(sqDoyCount,:) = counts;
            end
        end
    end
    % sunsetHist has to be 0 or 1
    odbaMeanSunset = [];
    % % odbaMaxMeanSunset = [];
    odbaIdsSunset = [];
    odbaMeanSunrise = [];
    % % odbaMaxMeanSunrise = [];
    odbaIdsSunrise = [];
    for iBin = 1:size(sunsetHist,2)
        useIds_odbaSunset = find(sunsetHist(:,iBin));
        odbaMeanSunset(iBin) = mean(odbaHist(useIds_odbaSunset));
        % %     odbaMaxMeanSunset(iBin) = mean(odbaMaxHist(useIds_odbaSunset));
        odbaIdsSunset(iBin) = numel(useIds_odbaSunset);
        
        useIds_odbaSunrise = find(sunriseHist(:,iBin));
        odbaMeanSunrise(iBin) = mean(odbaHist(useIds_odbaSunset));
        % %     odbaMaxMeanSunrise(iBin) = mean(odbaMaxHist(useIds_odbaSunset));
        odbaIdsSunrise(iBin) = numel(useIds_odbaSunrise);
    end
    
    % f = fit(odbaMeanSunset',odbaIdsSunset','poly1');
    % figure; plot(f,odbaMeanSunset',odbaIdsSunset');
    %
    % f = fit(odbaMeanSunrise',odbaIdsSunrise','poly1');
    % figure; plot(f,odbaMeanSunrise',odbaIdsSunrise');
    
    %% plot
    thisMeanSunset = [];
    thisMaxMeanSunset = [];
    thisIdsSunset = [];
    thisMeanSunrise = [];
    thisMaxMeanSunrise = [];
    thisIdsSunrise = [];
    usePoints = 30; % half window
    for ii = usePoints:numel(odbaIdsSunrise)-usePoints
        nR = ii-usePoints+1:ii+usePoints;
        thisMeanSunset(ii,:) = odbaMeanSunset(nR);
        % %     thisMaxMeanSunset(ii,:) = odbaMaxMeanSunset(nR);
        thisIdsSunset(ii,:) = odbaIdsSunset(nR);
        thisMeanSunrise(ii,:) = odbaMeanSunrise(nR);
        % %     thisMaxMeanSunrise(ii,:) = odbaMaxMeanSunrise(nR);
        thisIdsSunrise(ii,:) = odbaIdsSunrise(nR);
    end
    
    % close all
    
    % ff(1200,600);
    % subplot(211);
    % rs = [];
    % ps = [];
    % for ii = 1:size(thisMeanSunset,1)
    %     [r,p] = corr(thisMeanSunset(ii,:)',thisIdsSunset(ii,:)');
    %     rs(ii) = r;
    %     ps(ii) = p;
    % end
    % plot(t,rs);
    % yyaxis right;
    % plot(t,ps);
    % ylim([0 0.1]);
    % xlim([min(t),max(t)]);
    % xticklabels(xticks/3600);
    % title('sunset');
    
    subplot(4,1,iSeason);
    rs = [];
    ps = [];
    for ii = 1:size(thisMeanSunrise,1)
        [r,p] = corr(thisMeanSunrise(ii,:)',thisIdsSunrise(ii,:)');
        rs(ii) = r;
        ps(ii) = p;
    end
    t = binEdges((1:numel(rs))+usePoints-1);
    plot(t,rs);
    yyaxis right;
    plot(t,ps);
    ylim([0 0.1]);
    xlim([min(t),max(t)]);
    xticklabels(xticks/3600);
    title(seasonLabel{iSeason});
end


% this is a good view of whats going on with the data
% % fs = 12;
% % close all
% % ff(1200,800);
% % for ii = 1:2
% %     subplot(2,1,ii);
% %     if ii == 1
% %         titleLabel = 'Awake-sunset';
% %         binEdges = linspace(0,nHours*3600,nBins);
% %         counts = histcounts(awake_sunset,binEdges);
% %     else
% %         titleLabel = 'Awake-sunrise';
% %         binEdges = -fliplr(linspace(0,nHours*3600,nBins));
% %         counts = histcounts(awake_sunrise,binEdges);
% %     end
% %     t = binEdges(1:end-1);
% %     plot(t,counts,'k-');
% %     title(titleLabel);
% %     xlim([min(t) max(t)]);
% %     xticklabels(compose('%1.2f',xticks/3600));
% %     xlabel('rel. hours');
% %     if ii == 1
% %         text(min(xlim),1,'\leftarrow sunset','fontsize',fs);
% %     else
% %         text(max(ylim),1,'\rightarrow sunrise','fontsize',fs);
% %     end
% % end