sqkey = readtable('sqkey');
filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
nBins = 96;
nInt = 1440/96;
squirrelId = 0;
sIds = round(linspace(1,366,5));
seasonDoys = circshift(1:366,60);
if do
    cosStats = [];
    seasonMap = [];
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq}) % && ~any(ismember(sqkey.year(iSq),[2014,2019])) % ~(strcmp(sqkey.source{iSq},'ES') &&
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            T = detect_sleepWake(T);
            squirrelId = squirrelId + 1;
            dtdoys = day(T.datetime,'dayofyear');
            undoys = unique(dtdoys);
            for iSeason = 1:4
                if ismember(undoys(round(numel(undoys)/2)),seasonDoys(sIds(iSeason):sIds(iSeason+1)))
                    seasonMap(squirrelId) = iSeason;
                end
            end
            actograms = [];
            Ts = smoothdata(movsum(T.awake,60),'gaussian',1);
            for iDoy = 1:numel(undoys)
                theseDoys = find(dtdoys == undoys(iDoy));
                if numel(theseDoys) == 1440
                    sunrise = Tss.sunrise(Tss_doys == undoys(iDoy));
                    afterSunrise = Tss.sunset(Tss_doys == undoys(iDoy));
                    closestId = closest(secDay(T.datetime(theseDoys)),secDay(sunrise)); % center on sunrise
                    theseDoys = theseDoys(closestId):theseDoys(closestId) + 1440-1;
                    if min(theseDoys) > 1 && max(theseDoys) < numel(T.datetime)
                        actograms = [actograms;Ts(theseDoys)];
                    end
                end
            end
% % %     %         y = smoothdata(movsum(T.awake,60),'gaussian',1)';
            if isempty(actograms)
                squirrelId = squirrelId - 1;
            else
                y = normalize((actograms),'range')'; % 0-1 interval
                t = linspace(0,2*pi,numel(y));%1:numel(y);
                w = numel(y)/1440;
                alpha = .05;
                close all;

        %         t = linspace(0,2*pi,1000);
        %         y = sin(t);
        %         w = 1000;
                [M,Amp,phi] = cosinor(t,y,w,alpha); % M,Amp,phi
                cosStats(squirrelId,:) = [M,Amp,phi];
            end
        end
    end
    do = false;
end
%%
clc
sTitles = {'winter','spring','summer','fall','all'};
for iSeason = 1:5
    if iSeason == 5
        useIds = 1:size(cosStats,1);
    else
        useIds = seasonMap == iSeason;
    end
    disp(sTitles{iSeason});
    fprintf('Mesor: %1.2f + %1.2f \nAmp: %1.2f + %1.2f \nAcrophase: %1.2f + %1.2f (%1.2f hr + %1.2f hr)\n\n',...
        mean(cosStats(useIds,1)),std(cosStats(useIds,1)),...
        mean(cosStats(useIds,2)),std(cosStats(useIds,2)),...
        mean(cosStats(useIds,3)),std(cosStats(useIds,3)),...
        24*mean(cosStats(useIds,3))/(2*pi),24*std(cosStats(useIds,3))/(2*pi));
end