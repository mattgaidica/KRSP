function T_a = find_dayNight(loadfile,doFig,doWrite)
nSmooth = 360;
colors = lines(3);
T_a = table;

load(loadfile);

odba_gauss = smoothdata(smooth(T.odba,nSmooth/8),'gaussian',nSmooth);
odba_z = (odba_gauss - mean(odba_gauss)) / std(odba_gauss);
odba_bin = zeros(size(odba_z));
zThresh = 0;
odba_bin(odba_z > zThresh) = 1; % set thresh here

all_awake_ids = find(diff(odba_bin) > 0);
all_asleep_ids = find(diff(odba_bin) < 0);

allDays = unique(datetime(year(T.datetime),month(T.datetime),day(T.datetime)));
sunrise_id = NaN(numel(allDays),1);
sunset_id = NaN(numel(allDays),1);
awake_id = NaN(numel(allDays),1);
asleep_id = NaN(numel(allDays),1);
for iDay = 1:numel(allDays)
    [sunrise,sunset,day_length] = sunriseSunset(allDays(iDay));
    % find the closest timestamp in the decimated data
    tryId = find(T.datetime - sunrise > 0,1,'first');
    T_a.sunrise(iDay) = sunrise;
    if ~isempty(tryId)
        sunrise_id(iDay) = tryId;
        % find the closest awake/asleep identifier
        [~,id] = closest(all_awake_ids,sunrise_id(iDay));
        T_a.awake(iDay) = T.datetime(id);
        T_a.awake_sunrise(iDay) = seconds(T_a.awake(iDay) - sunrise);
        awake_id(iDay) = id;
    end
    
    tryId = find(T.datetime - sunset > 0,1,'first');
    T_a.sunset(iDay) = sunset;
    if ~isempty(tryId)
        sunset_id(iDay) = tryId;
        [~,id] = closest(all_asleep_ids,sunset_id(iDay));
        T_a.asleep(iDay) = T.datetime(id);
        T_a.asleep_sunset(iDay) = seconds(T_a.asleep(iDay) - sunset);
        asleep_id(iDay) = id;
    end
    T_a.day_length(iDay) = day_length;
end

if doWrite
    save(strrep(loadfile,'.csv.mat','.csv.meta.mat'),'T_a','nSmooth','zThresh');
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
    for iDay = 1:numel(allDays)
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
    xlabel_ids = sort([sunrise_id;sunset_id]);
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
    for iDay = 1:numel(allDays)
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