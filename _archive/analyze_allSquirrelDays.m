load('sqkey.mat');
dataPath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
warning('off','all');

useDoy = false;
colors = lines(5);

if do
    uniSq = 0;
    Td = table;
    useSquirrels = [];
    for ii = 1:size(sqkey,1)
        if ~isempty(sqkey.filenames{ii}) && isfile(fullfile(dataPath,sqkey.filenames{ii}))...
                && any(strcmp(sqkey.Sex{ii},{'M','F'}))
            fn = sqkey.filenames{ii};
            fn_meta = strrep(fn,'.mat','_meta.mat');
            % manual exclusion
            if any(strcmp(fn_meta,{'E11_Mar10_2016__20160325_meta.mat'}))
                continue;
            end
            uniSq = uniSq + 1;
            useSquirrels(uniSq) = ii;
            disp(fn_meta);
            load(fullfile(dataPath,fn_meta));
            for iDay = 1:size(T_a,1)
                if useDoy
                    dn = day(T_a.sunrise(iDay),'dayofyear');
                else
                    dn = datetime(year(T_a.sunrise(iDay)),month(T_a.sunrise(iDay)),day(T_a.sunrise(iDay)));
                end
                [~,~,~,H,MN,S] = datevec(T_a.awake(iDay));
                ds_awake = H*3600+MN*60+S;
                [~,~,~,H,MN,S] = datevec(T_a.asleep(iDay));
                ds_asleep = H*3600+MN*60+S;
                if isempty(Td)
                    useRow = 1;
                    Td.day = dn; % init
                    Td.awakes(useRow) = {ds_awake};
                    Td.asleeps(useRow) = {ds_asleep};
                    Td.awake_odba(useRow) = {T_a.awake_odba(iDay)};
                    Td.asleep_odba(useRow) = {T_a.asleep_odba(iDay)};
                    norm_odba = T_a.awake_odba(iDay) / (ds_asleep - ds_awake);
                    Td.norm_odba(useRow) = {norm_odba};
                else
                    useRow = find(Td.day == dn,1);
                    if isempty(useRow)
                        useRow = numel(Td.day)+1;
                        Td.day(useRow) = dn;
                        Td.awakes(useRow) = {ds_awake};
                        Td.asleeps(useRow) = {ds_asleep};
                        Td.awake_odba(useRow) = {T_a.awake_odba(iDay)};
                        Td.asleep_odba(useRow) = {T_a.asleep_odba(iDay)};
                        norm_odba = T_a.awake_odba(iDay) / (ds_asleep - ds_awake);
                        Td.norm_odba(useRow) = {norm_odba};
                    else
                        Td.awakes(useRow) = {[Td.awakes{useRow} ds_awake]};
                        Td.asleeps(useRow) = {[Td.asleeps{useRow} ds_asleep]};
                        Td.awake_odba(useRow) = {[Td.awake_odba{useRow} T_a.awake_odba(iDay)]};
                        Td.asleep_odba(useRow) = {[Td.asleep_odba{useRow} T_a.asleep_odba(iDay)]};
                        norm_odba = T_a.awake_odba(iDay) / (ds_asleep - ds_awake);
                        Td.norm_odba(useRow) = {[Td.norm_odba{useRow} norm_odba]};
                    end
                end
            end
        end
    end
    do = false;
end

% if ~exist('ss_x')
    ss_x = table;
    if useDoy
        for iDay = 1:366%numel(Td.day)
            curDay = datetime('1-Jan-2016') + iDay - 1; % includes leap year
            ss_x.date(iDay) = iDay;
            [sunrise,sunset,day_length] = sunriseSunset(curDay);
            [~,~,~,H,MN,S] = datevec(sunrise);
            ss_x.sunrise(iDay) = H*3600+MN*60+S;
            [~,~,~,H,MN,S] = datevec(sunset);
            ss_x.sunset(iDay) = H*3600+MN*60+S;
            TdId = find(Td.day == iDay);
            if ~isempty(TdId)
                ss_x.norm_odba(iDay) = nanmean(Td.norm_odba{TdId});
            else
                ss_x.norm_odba(iDay) = NaN;
            end
        end
    else
        for iDay = 1:days(max(Td.day)-min(Td.day))
            curDay = min(Td.day) + days(iDay-1);
            ss_x.date(iDay) = curDay;
            [sunrise,sunset,day_length] = sunriseSunset(curDay);
            [~,~,~,H,MN,S] = datevec(sunrise);
            ss_x.sunrise(iDay) = H*3600+MN*60+S;
            [~,~,~,H,MN,S] = datevec(sunset);
            ss_x.sunset(iDay) = H*3600+MN*60+S;
            TdId = find(Td.day == curDay);
            if ~isempty(TdId)
                ss_x.norm_odba(iDay) = nanmean(Td.norm_odba{TdId});
            else
                ss_x.norm_odba(iDay) = NaN;
            end
        end
%     end
end

lns = [];
op = 0.1;
close all
ff(1200,600);
% for iDay = 1:numel(Td.day)
%     plot([Td.day(iDay) Td.day(iDay)],[mean(Td.awakes{iDay})+std(Td.awakes{iDay}),...
%         mean(Td.awakes{iDay})-std(Td.awakes{iDay})],'-','color',[0 0 0 op]);
%     plot([Td.day(iDay) Td.day(iDay)],[mean(Td.asleeps{iDay})+std(Td.asleeps{iDay}),...
%         mean(Td.asleeps{iDay})-std(Td.asleeps{iDay})],'-','color',[0 0 0 op]);
%     hold on;
% end
lns(1) = plot(Td.day,cellfun(@nanmean,Td.awakes),'.','color','r');
hold on;
lns(2) = plot(Td.day,cellfun(@nanmean,Td.asleeps),'x','color','r');
ylim([0 24*3600]);
yticks([0:3*3600:24*3600]);
yticklabels(compose('%1.2f',yticks/3600));
ylabel('hour of day');

if useDoy
    lns(5) = plot(ss_x.date,ss_x.sunrise,'-','color',[colors(3,:) 0.25],'linewidth',3);
    lns(6) = plot(ss_x.date,ss_x.sunset,'-','color',[0 0 0 0.25],'linewidth',3);
    xlim([1 366]);
    xlabel('day of year');
else
    lns(5) = plot(ss_x.date,ss_x.sunrise,'-','color',[colors(3,:) 0.25],'linewidth',3);
    lns(6) = plot(ss_x.date,ss_x.sunset,'-','color',[0 0 0 0.25],'linewidth',3);
    xlim([min(Td.day) max(Td.day)]);
end
yyaxis right;
lns(3) = bar(Td.day,cellfun(@nanmean,Td.awake_odba),'edgecolor','none','facecolor','k');
hold on;
lns(4) = bar(Td.day,-cellfun(@nanmean,Td.asleep_odba),'edgecolor','none','facecolor',colors(3,:));
lns(7) = plot(ss_x.date,smoothdata(ss_x.norm_odba,'gaussian',25)*1e5,'-','color',colors(5,:));
ylim([-800,800]);
ylabel('sum odba');
set(gca,'ycolor','k');

yyaxis left;
if useDoy
    for iMonth = 1:12
        dt = datetime(2016,iMonth,1);
        doy = day(dt,'dayofyear');
        plot([doy doy],[min(ylim) max(ylim)],'k:');
        text(doy,min(ylim)+3600,['\leftarrow',datestr(dt,'mmm')]);
    end
else
    all_years = unique(year(ss_x.date));
    for iYear = 1:numel(all_years)
        janId = find(month(ss_x.date) == 1 & year(ss_x.date) == all_years(iYear),1,'first');
        if ~isempty(janId)
            plot([ss_x.date(janId) ss_x.date(janId)],[min(ylim) max(ylim)],'k:');
        end
    end
end

legend(lns,{'awake','asleep','odba day','odba night','sunrise','sunset','norm-scaled odba'},...
    'location','eastoutside');
set(gca,'fontsize',20);
warning('on','all');