function T_a = find_dayNight(loadfile,doFig,doWrite)
warning('off','all');
nSmooth = 360;
colors = lines(3);
T_a = table;

load(loadfile);

odba = sqrt(smooth(T.odba,nSmooth/8));
odba_gauss = smoothdata(odba,'gaussian',nSmooth);
odba_z = (odba_gauss - mean(odba_gauss)) / std(odba_gauss);
odba_bin = zeros(size(odba_z));
zThresh = 0;
odba_bin(odba_z > zThresh) = 1; % set thresh here

all_awake_ids = find(diff(odba_bin) > 0);
all_asleep_ids = find(diff(odba_bin) < 0);

allDays = unique(datetime(year(T.datetime),month(T.datetime),day(T.datetime)));
dayCount = 0;
sunrise_id = [];
sunset_id = [];
awake_id = [];
asleep_id = [];
for iDay = 1:numel(allDays)
    [sunrise,sunset,day_length] = sunriseSunset(allDays(iDay));
    % find the closest timestamp in the decimated data
    tryId_sunrise = find(T.datetime - sunrise > 0,1,'first');
    tryId_sunset = find(T.datetime - sunset > 0,1,'first');
    if ~isempty(tryId_sunrise) && ~isempty(tryId_sunset)
        dayCount = dayCount + 1;
        T_a.sunrise(dayCount) = sunrise;
        if ~isempty(tryId_sunrise)
            sunrise_id(dayCount) = tryId_sunrise;
            % find the closest awake/asleep identifier
            [~,id] = closest(all_awake_ids,sunrise_id(dayCount));
            T_a.awake(dayCount) = T.datetime(id);
            T_a.awake_sunrise(dayCount) = seconds(T_a.awake(dayCount) - sunrise);
            awake_id(dayCount) = id;
        end
        
        T_a.sunset(dayCount) = sunset;
        if ~isempty(tryId_sunset)
            sunset_id(dayCount) = tryId_sunset;
            [~,id] = closest(all_asleep_ids,sunset_id(dayCount));
            T_a.asleep(dayCount) = T.datetime(id);
            T_a.asleep_sunset(dayCount) = seconds(T_a.asleep(dayCount) - sunset);
            asleep_id(dayCount) = id;
        end
        if ~isempty(awake_id(dayCount)) && ~isempty(asleep_id(dayCount))
            T_a.awake_odba(dayCount) = sum(T.odba(awake_id(dayCount):asleep_id(dayCount)));
            if dayCount > 1
                T_a.asleep_odba(dayCount-1) = sum(T.odba(asleep_id(dayCount-1):awake_id(dayCount)));
            end
            if dayCount == numel(allDays)
                T_a.asleep_odba(dayCount) = NaN; % always end with partial data?
            end
        end
        T_a.day_length(dayCount) = day_length;
    end
end
T_a.awake_odba(T_a.awake_odba==0) = NaN;
T_a.asleep_odba(T_a.asleep_odba==0) = NaN;

if doWrite
    save(strrep(loadfile,'.mat','_meta.mat'),'T_a','nSmooth','zThresh');
end

if doFig
    op = 0.15;
    fs = 10;
    ly = 1.15;
    lns = [];
    close all
    h = ff(1900,800);
    
    yyaxis right;
    hold on;
    for iDay = 1:dayCount
        if ~isempty(sunrise_id(iDay))
            lns(5) = plot([sunrise_id(iDay) sunrise_id(iDay)],ylim,'-','color',[colors(3,:)],'linewidth',4);
        end
        if ~isempty(sunset_id(iDay))
            lns(6) = plot([sunset_id(iDay) sunset_id(iDay)],ylim,'-','color',[0,0,0],'linewidth',4);
        end
    end
    lns(2) = plot(odba_z,'k-','lineWidth',2);
    hold on;
    ylabel('z-odba');
    set(gca,'ycolor','k')
    xlabel_ids = sort([sunrise_id sunset_id]);
    xticks(xlabel_ids);
    xticklabels(datestr(T.datetime(xlabel_ids)));
    xtickangle(20);
    xlabel('sunrise/sunset');
    xlim([1 numel(odba_z)]);
    set(gca,'fontsize',fs);
    lns(3) = plot(odba_bin,'--r');
    lns(4) = plot(xlim,[zThresh,zThresh],':r','linewidth',2);
    
    yyaxis left;
    lns(1) = plot(T.odba,'color',[colors(1,:) op]);
    ylabel('raw odba');
    set(gca,'ycolor',colors(1,:))
    set(gca,'fontsize',fs);
    
    yyaxis right;
    for iDay = 1:dayCount
        if ~isempty(awake_id(iDay))
            text(awake_id(iDay),ly+0.05,{['#',num2str(iDay)],'\uparrow'},'fontsize',fs-2,'horizontalalign','center','color','r');
        end
        if ~isempty(asleep_id(iDay))
            text(asleep_id(iDay),ly,{['#',num2str(iDay)],'\downarrow'},'fontsize',fs-2,'horizontalalign','center','color','r');
        end
    end
    legend(lns,{'raw odba','gaussian z-odba','binary classifier','z-thresh','sunrise','sunset'},'location','eastoutside','fontsize',fs);
    saveas(h,[loadfile,'.png']);
    close(h);
end
warning('on','all');