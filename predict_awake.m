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
            if ~isValidT(T,false) % no temp
                disp('skipping');
                continue;
            end
        end
        disp(sqkey.filename{iSq});
        T = detect_sleepWake(T);
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

% close all
% % % % useMonths = [1,3;4,6;7,9;10,12];
% % % % useDoys = [1,30;70,110;170,190;260,280];
sIds = round(linspace(1,366,5));
seasonDoys = circshift(1:366,60);
% % % % useDoys = {[350:366 1:20],70:110,160:190,235:295};
useDoys = {seasonDoys(sIds(1):sIds(2)),seasonDoys(sIds(2):sIds(3)),...
    seasonDoys(sIds(3):sIds(4)),seasonDoys(sIds(4):sIds(5))};
monthNames = {'Winter','Spring','Summer','Fall'};
typeNames = {'asleep','awake'};
ff(750,600,2);
rows = size(useDoys,2);
cols = 2;
nS = 1;
for iType = 0:1
    for iiDoy = 1:size(useDoys,2)
        subplot(rows,cols,prc(cols,[iiDoy,iType+1]));
        binEdges_x = linspace(0,86400,24*4+1);
        binEdges_y = logspace(2.9,5,48);
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
        imagesc(imgaussfilt(trans_mat' ./ trans_norm',nS));
        hold on;
        set(gca,'ydir','normal');
        xticks(linspace(min(xlim),max(xlim),8));
        xticklabels(compose('%1.0f',linspace(0,23,numel(xticks))));
        if iiDoy == size(useDoys,2)
            xlabel('hour of day');
        end
        caxis([0 0.15]);
        yticklabels(compose('%1.1f',binEdges_y(yticks)/60/60));
        if iType == 0
            ylabel('hours');
        end
        c = colorbar;
        c.Label.String = ['p(',typeNames{iType+1},')'];
        c.Label.FontSize = 14;
        title(monthNames{iiDoy});
        colormap(magma);
        
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