loadPath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
files = dir(fullfile(loadPath,'*.mat'));
strFilts = {'Nest','Out'};
doDebug = false;
op = 0.2;
colors = lines(5);
if do
    Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
    Tss_doys = day(Tss.sunrise,'dayofyear');
    trans_light = [];
    trans_dark = [];
    trans_all = [];
    trans_secs = [];
    trans_months = [];
    trans_type = [];
    for iFile = 1:numel(files)
        disp(files(iFile).name);
        load(fullfile(loadPath,files(iFile).name));
        filtDt = find(seconds(diff(Tstat.datetime)) < 60);
        filtDt = [filtDt;filtDt+1];
        Tstat(filtDt,:) = [];
        
        if doDebug
            ff(1200,600);
            Tnest = strcmp(T.nest,'Nest');
            plot(T.datetime,T.odba,'k-');
            ylabel('ODBA');
            hold on;
            plot(T.datetime,Tnest,'r-');
            yyaxis right;
            plot(T.datetime,T.temp,'m-');
            set(gca,'ycolor','m')
            ylabel('temp (C)');
        end
        
        for iT = 1:size(Tstat,1)-1 % last entry has no transition to reference
            doy = day(Tstat.datetime(iT),'dayofyear');
            % does this entry occur in light or dark?
            sunrise = Tss.sunrise(Tss_doys == doy);
            sunset = Tss.sunset(Tss_doys == doy);
            transTime = seconds(Tstat.datetime(iT+1) - Tstat.datetime(iT));
            trans_all = [trans_all transTime];
            trans_secs = [trans_secs secDay(Tstat.datetime(iT))];
            trans_months = [trans_months month(Tstat.datetime(iT))];
            trans_type = [trans_type strcmp(Tstat.nest(iT),'Nest')];
            if doDebug
                yyaxis left;
                plot(Tstat.datetime(iT),1,'g*');
                text(Tstat.datetime(iT),(mod(iT,6)/3)+1,[num2str(transTime),'s'],'color','r');
                tdt = Tstat.datetime(iT);
                sunrise_x = datetime(year(tdt),month(tdt),day(tdt),hour(sunrise),minute(sunrise),0);
                sunset_x = datetime(year(tdt),month(tdt),day(tdt),hour(sunset),minute(sunset),0);
                plot([sunrise_x,sunrise_x],ylim,'-','linewidth',4,'color',[colors(3,:),op]);
                plot([sunset_x,sunset_x],ylim,'-','linewidth',4,'color',[0 0 0 op]);
            end
            if secDay(Tstat.datetime(iT)) >= secDay(sunrise) && secDay(Tstat.datetime(iT)) < secDay(sunset)
                trans_light = [trans_light;transTime];
            else
                trans_dark = [trans_dark;transTime];
            end
        end
    end
    do = false;
end

Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
Tss_months = month(Tss.sunrise);

close all
useMonths = [1,3;4,6;7,9;10,12];
% useMonths = [1:12;1:12]';
monthNames = {'Winter','Spring','Summer','Fall'};
typeNames = {'in','out'};
ff(1400,900);
rows = size(useMonths,1);
cols = 2;
for iType = 0:1
    for iMonth = 1:size(useMonths,1)
        subplot(rows,cols,prc(cols,[iMonth,iType+1]));
        binEdges_x = linspace(0,86400,24*4+1);
        %         binEdges_y = linspace(0,12*60*60,40+1); % H*M*S
        binEdges_y = logspace(1,5,80);
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
        c.Label.String = ['p(transition ',typeNames{iType+1},')'];
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
        plot([sunrise_std_x,sunset_std_x],[3,3],'linewidth',8,'color',[colors(3,:) 0.5]);
        plot([sunrise_x,sunset_x],[3,3],'linewidth',8,'color',[colors(3,:) 0.8]);
        
        set(gca,'fontsize',14);
        set(gca,'TitleFontSizeMultiplier',1.25);
    end
end

binEdges = linspace(0,8*60*60,20); % 10min to 5h
% close all
ff(1200,600,2);

subplot(121);
histogram(trans_light,binEdges);
% xticklabels(compose('%1.2f',xticks/3600));
xticklabels(compose('%1.0f',xticks));
xlabel('s');
title('light');
% ylim([0 3200]);

subplot(122);
histogram(trans_dark,binEdges);
% xticklabels(compose('%1.2f',xticks/3600));
xticklabels(compose('%1.0f',xticks));
xlabel('s');
title('dark');
% ylim([0 3200]);