if do
    sqkey = readtable('sqkey');
    filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
    Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
    Tss_doys = day(Tss.sunrise,'dayofyear');
    sq_odba = [];
    sq_awake = [];
    sq_ids = [];
    sq_doys = [];
    sq_dayLength = [];
    iRow = 0;
    squirrelId = 0;
    tWindow = 4*60; % 4 hours (min sunrise ~= 4.5hrs)
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq}) && ~any(ismember(sqkey.year(iSq),[2014,2019])) % ~(strcmp(sqkey.source{iSq},'ES') && 
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            T = detect_sleepWake(T,2);
            if isValidT(T,false)
                dtdoys = day(T.datetime,'dayofyear');
                undoys = unique(dtdoys);
                squirrelId = squirrelId + 1;
                for iDoy = 1:numel(undoys)
                    if any(ismember(271:279,undoys(iDoy))) && strcmp(sqkey.source{iSq},'BD')
                        continue;
                    end
                    theseDoys = find(dtdoys == undoys(iDoy));
                    if numel(theseDoys) == 1440 % require full day for now
                        iRow = iRow + 1;
                        sunrise = Tss.sunrise(Tss_doys == undoys(iDoy));
                        sunset = Tss.sunset(Tss_doys == undoys(iDoy));
                        sq_ids(iRow) = squirrelId;
                        sq_doys(iRow) = undoys(iDoy);
                        sq_odba(iRow,:) = T.odba(theseDoys);
                        sq_dayLength(iRow) = Tss.day_length(Tss_doys == undoys(iDoy));
                        closestId = closest(secDay(T.datetime(theseDoys)),secDay(sunrise)); % center on sunrise
                        sq_awake(iRow,:) = T.awake(closestId-tWindow:closestId+tWindow-1);
                    end
                end
            end
        end
    end
    do = false;
end

% % norm_method = {'zscore','norm','scale','range'};
% % norm_caxis = [-0.25 1.1;0 0.045;0 2.05;0 0.175];
% % rows = 4;
% % cols = 2;
% % mean_odba = [];
% % mean_norm_odba = [];
% % close all;
% % ff(1200,1200);
% % iSubplot = 0;
% % for iNorm = 1:4
% %     for iDoy = 1:366
% %         useIds = find(sq_doys == iDoy);
% %         mean_odba(iDoy,:) = mean(sq_odba(useIds,:));
% %         if ~isempty(useIds)
% %             normvals = [];
% %             for ii = 1:numel(useIds)
% %                 thisOdba = sq_odba(useIds(ii),:);
% %                 normvals(ii,:) = normalize(thisOdba,norm_method{iNorm});
% %             end
% %             mean_norm_odba(iDoy,:) = mean(normvals);
% %         end
% %     end
% %     iSubplot = iSubplot + 1;
% %     subplot(rows,cols,iSubplot);
% %     imagesc(mean_odba'); set(gca,'ydir','normal');
% %     xlabel('doy');ylabel('min. into day');
% %     colorbar;colormap(magma);caxis([0 1.3]);
% %     title('ODBA');
% %     
% %     iSubplot = iSubplot + 1;
% %     subplot(rows,cols,iSubplot);
% %     imagesc(mean_norm_odba'); set(gca,'ydir','normal');
% %     xlabel('doy');ylabel('min. into day');
% %     colorbar;colormap(magma);caxis(norm_caxis(iNorm,:));
% %     title(['norm ODBA by ',norm_method{iNorm}]);
% % end


mean_odba = [];
mean_norm_odba = zeros(366,1440);
max_odba = [];
daylength_odba = [];
mean_awake = zeros(366,480);
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
        mean_norm_odba(iDoy,:) = nanmean(normvals);
    end
%     daylength_odba(iDoy) = mean_odba(iDoy) ./ mean(sq_dayLength(useIds));
end
% close all
t = linspace(-tWindow,tWindow,size(sq_awake,2));
ff(800,1100);
subplot(311);
imagesc(mean_odba'); set(gca,'ydir','normal');
xlabel('doy');ylabel('min. into day');
colorbar;colormap(magma);caxis([0 1]);
title('ODBA');

subplot(312);
imagesc(mean_norm_odba'); set(gca,'ydir','normal');
xlabel('doy');ylabel('min. into day');
colorbar;colormap(magma);%caxis([-.25 1.5]);
title('norm ODBA');

subplot(313);
imagesc(1:366,t,mean_awake'); set(gca,'ydir','normal');
xlabel('doy');ylabel('min. from sunrise');
colorbar;colormap(magma);
title('mean awake');

% plot(mean_awake);

t = linspace(-tWindow,tWindow,size(sq_awake,2));
mean_odba = [];
z_mean_odba = [];
nSurr = 3;
seasonDays = round(linspace(1,366,5));
iSeason = 4;
for iBin = 1:size(sq_awake,2)
    awakeIds = find(sq_awake(:,iBin)==1 &...
        (sq_doys' >= seasonDays(iSeason) & sq_doys' <= seasonDays(iSeason+1))...
    );
    mean_odba(iBin) = mean(sq_odba(awakeIds));
    for iSurr = 1:nSurr
        surr_awakeIds = randi([1,size(sq_awake,1)],numel(awakeIds),1);
        surr_mean_odba(iSurr) = mean(sq_odba(surr_awakeIds));
    end
    z_mean_odba(iBin) = (mean_odba(iBin) - mean(surr_mean_odba)) / std(surr_mean_odba);
end

% close all
% ff(900,700);
subplot(311);
plot(t,sum(sq_awake)); hold on;
subplot(312);
plot(t,mean_odba);
hold on;
plot(xlim,[mean(sq_odba) mean(sq_odba)],'r-');
subplot(313);
plot(t,z_mean_odba); hold on