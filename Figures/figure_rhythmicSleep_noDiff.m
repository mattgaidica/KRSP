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
    tWindow = 4*60; % 4 hours (min sunrise ~= 4.5hrs)
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq})
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            T = detect_sleepWake(T,2);
            dtdoys = day(T.datetime,'dayofyear');
            undoys = unique(dtdoys);
            squirrelId = squirrelId + 1;
            for iDoy = 1:numel(undoys)
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
    do = false;
end

mean_odba = [];
max_odba = [];
daylength_odba = [];
mean_awake = [];
for iDoy = 1:366
    useIds = find(sq_doys == iDoy);
%     mean_odba(iDoy) = sq_odba(useIds,:);
    if ~isempty(useIds)
        mean_awake(iDoy,:) = median(sq_awake(useIds,:),1);
        maxvals = [];
        for ii = 1:numel(useIds)
            thisOdba = sq_odba(useIds(ii),:);
            v = sort(thisOdba);
            maxvals(ii) = std(thisOdba);
        end
        max_odba(iDoy) = mean(maxvals);
    end
%     daylength_odba(iDoy) = mean_odba(iDoy) ./ mean(sq_dayLength(useIds));
end

t = linspace(-tWindow,tWindow,size(sq_awake,2));
ff(800,800);
subplot(211);
plot(mean_odba);
subplot(212);
imagesc(1:366,t,mean_awake'); set(gca,'ydir','normal');
xlabel('doy');ylabel('min. from sunrise');
colorbar;colormap(magma);
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