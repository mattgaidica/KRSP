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
    trans_type = [];
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq})
            load(fullfile(loadPath,sqkey.filename{iSq}));
            if ~isValidT(T,false) % no temp
                continue;
            end
        end
        disp(sqkey.filename{iSq});
        T = detect_sleepWake(T,2);
        Tawake = make_Tawake(T);
        
        for iT = 1:size(Tawake,1)-1 % last entry has no transition to reference
            transTime = seconds(Tawake.datetime(iT+1) - Tawake.datetime(iT));
            trans_all = [trans_all transTime];
            trans_secs = [trans_secs secDay(Tawake.datetime(iT))];
            trans_months = [trans_months month(Tawake.datetime(iT))];
            trans_type = [trans_type Tawake.awake(iT)];
        end
    end
    do = false;
end

Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
Tss_months = month(Tss.sunrise);

close all
useMonths = [1,3;4,6;7,9;10,12];
monthNames = {'Winter','Spring','Summer','Fall'};
typeNames = {'awake','asleep'};
ff(900,900);
rows = size(useMonths,1);
cols = 2;
for iType = 0:1
    for iMonth = 1:size(useMonths,1)
        subplot(rows,cols,prc(cols,[iMonth,iType+1]));
        binEdges_x = linspace(0,86400,24*2+1);
        binEdges_y = logspace(2.8,4.8,24);
        trans_mat = zeros(numel(binEdges_x)-1,numel(binEdges_y)-1);
        trans_norm = zeros(numel(binEdges_x)-1,1);
        for iBin = 1:numel(binEdges_x)-1
            useIds = find(trans_secs >= binEdges_x(iBin) & trans_secs < binEdges_x(iBin+1)...
                & trans_months >= useMonths(iMonth,1) & trans_months <= useMonths(iMonth,2)...
                & trans_type == iType);
            theseVals = trans_all(useIds);
            n = histcounts(theseVals,binEdges_y);
            trans_mat(iBin,:) = n;
            trans_norm(iBin) = sum(n);
        end
        imagesc(trans_mat' ./ trans_norm');
        hold on;
        set(gca,'ydir','normal');
        xticklabels(compose('%1.0f',xticks/4));
        xlabel('hour of day');
        caxis([0 0.15]);
        yticklabels(compose('%1.0f',binEdges_y(yticks)/60));
        ylabel('minutes');
        c = colorbar;
        c.Label.String = ['p(',typeNames{iType+1},')'];
        c.Label.FontSize = 14;
        title(monthNames{iMonth});
        colormap(magma);
        
        tss_range = find(Tss_months >= useMonths(iMonth,1) & Tss_months <= useMonths(iMonth,2));
        theseSunrise = secDay(Tss.sunrise(tss_range));
        theseSunset = secDay(Tss.sunset(tss_range));
        xsc = linspace(0,86400,max(xlim)-min(xlim));
        sunrise_x = closest(xsc,mean(theseSunrise));
        sunrise_std_x = closest(xsc,mean(theseSunrise)-std(theseSunrise));
        sunset_x = closest(xsc,mean(theseSunset));
        sunset_std_x = closest(xsc,mean(theseSunset)+std(theseSunset));
        plot([sunrise_std_x,sunset_std_x],[1,1],'linewidth',8,'color',[colors(3,:) 0.5]);
        plot([sunrise_x,sunset_x],[1,1],'linewidth',8,'color',[colors(3,:) 0.8]);
        
        set(gca,'fontsize',14);
        set(gca,'TitleFontSizeMultiplier',1.25);
    end
end