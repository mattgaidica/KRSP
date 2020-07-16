if do
    sqkey = readtable('sqkey');
    filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
    Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
    Tss_doys = day(Tss.sunrise,'dayofyear');
    
    nBins = 24;
    binEdges = linspace(0,24,nBins+1) * 3600;
    
    doy_arr = [];
    sex_arr = [];
    odba_arr = zeros(1,nBins);
    asleep_arr = zeros(1,nBins);
    frac_out_nest = zeros(1,nBins);
    sqCount = 0;
    doyCount = 0;
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq})
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            T = detect_sleepWake(T,2);
            dtdoys = day(T.datetime,'dayofyear');
            undoys = unique(dtdoys);
            if numel(undoys) < 2
                continue;
            end
            sqCount = sqCount + 1;
            for iDoy = 1:numel(undoys)-1
                dts = T.datetime(dtdoys==undoys(iDoy) | dtdoys==undoys(iDoy+1));
                if (86400*2) - seconds(dts(end)-dts(1)) > 100 % not enough data
                    continue;
                end
                doyCount = doyCount + 1;
                theseOdba = T.odba(dtdoys==undoys(iDoy) | dtdoys==undoys(iDoy+1));
                theseNest = T.nest(dtdoys==undoys(iDoy) | dtdoys==undoys(iDoy+1));
                theseAsleep = ~T.awake(dtdoys==undoys(iDoy) | dtdoys==undoys(iDoy+1));
                theseSecs = secDay(dts);
                thisSunrise = secDay(Tss.sunrise(Tss_doys==undoys(iDoy)));
                [v,k] = sort(abs(theseSecs-thisSunrise));
                sunriseOdba = theseOdba(k(1):k(2));
                sunriseSleep = theseAsleep(k(1):k(2));
                sunriseSecs = theseSecs(k(1):k(2));
                sunriseSecs = sunriseSecs - sunriseSecs(1);
                addId = find(sunriseSecs<0,1,'first');
                sunriseSecs(addId:end) = sunriseSecs(addId:end) + sunriseSecs(addId-1) - sunriseSecs(addId);
                sunriseNest = theseNest(k(1):k(2));
                for iBin = 1:nBins
                    % note: mean of empty vector is NaN
                    tryIds = sunriseSecs >= binEdges(iBin) & sunriseSecs < binEdges(iBin+1);
                    if sum(tryIds) > 0
                        odba_arr(doyCount,iBin) = mean(sunriseOdba(tryIds));
                        asleep_arr(doyCount,iBin) = mean(sunriseSleep(tryIds));
                        frac_out_nest(doyCount,iBin) = sum(strcmp(sunriseNest(tryIds),'Out')) / sum(tryIds);
                    end
                end
                doy_arr(doyCount) = undoys(iDoy);
                sex_arr(doyCount) = strcmp(sqkey.sex{iSq},'M');
                mast_arr(doyCount) = ismember(sqkey.year(iSq),[2014,2019]);
            end
        end
    end
    do = false;
end

% % v = datevec(datenum(2016, ones(size(2016)), doy_arr));
% % mo_arr = v(:,2);

close all
op = 0.1;
sexLabels = {'F','M'};
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons.png',366);
ff(1400,800);
rows = 3;
cols = 3;
lineColors = lines(5);

for iSex = 0:1
    odba_compiled = zeros(366,nBins);
    asleep_compiled = zeros(366,nBins);
    frac_out_compiled = zeros(366,nBins);
    odba_compiled_std = zeros(366,nBins);
    for ii = 1:366
        tryIds = doy_arr==ii&sex_arr==iSex;
        if sum(tryIds) > 0
            odba_compiled(ii,:) = mean(odba_arr(tryIds,:));
            asleep_compiled(ii,:) = mean(asleep_arr(tryIds,:));
            odba_compiled_std(ii,:) = std(odba_arr(tryIds,:));
            frac_out_compiled(ii,:) = mean(frac_out_nest(tryIds,:));
        end
    end
    if iSex == 0
        odba_compiled_F = odba_compiled;
        asleep_compiled_F = asleep_compiled;
        odba_compiled_F_std = odba_compiled_std;
        frac_out_compiled_F = frac_out_compiled;
    else
        odba_compiled_M = odba_compiled;
        asleep_compiled_M = asleep_compiled;
        odba_compiled_M_std = odba_compiled_std;
        frac_out_compiled_M = frac_out_compiled;
    end
    
    % ODBA
    subplot(rows,cols,prc(cols,[iSex+1,1]));
    odba_cs = cumsum(odba_compiled')';
    useDoys = [];
    for ii = 1:366
        if sum(odba_cs(ii,:)) > 0% && numel(unique(diff(odba_cs(ii,:)))) > 10 % remove some flat lines
            useDoys = [useDoys;ii];
            plot(odba_cs(ii,:),'color',[colors(ii,:) op],'linewidth',2);
            hold on;
        end
    end
    drawnow;
    xlim([1 size(odba_cs,2)]);
    ylim([0 14]);
    yticks([0:2:max(ylim)]);
    xlabel('hours from sunrise');
    ylabel('cum sum ODBA');
    set(gca,'fontsize',14);
    title([sexLabels{iSex+1},' ','All (Out & In)']);
    cb = cbAside(gca,'doy','k',[1 366]);
    set(gca,'colormap',colors);
    grid on;
    
    % asleep
    subplot(rows,cols,prc(cols,[iSex+1,2]));
    asleep_cs = cumsum(asleep_compiled')';
    useDoys = [];
    for ii = 1:366
        if sum(asleep_cs(ii,:)) > 0% && numel(unique(diff(odba_cs(ii,:)))) > 10 % remove some flat lines
            useDoys = [useDoys;ii];
            plot(asleep_cs(ii,:),'color',[colors(ii,:) op],'linewidth',2);
            hold on;
        end
    end
    drawnow;
    xlim([1 size(asleep_cs,2)]);
    ylim([0 25]);
    yticks([0:2:max(ylim)]);
    xlabel('hours from sunrise');
    ylabel('cum sum asleep');
    set(gca,'fontsize',14);
    title([sexLabels{iSex+1},' ','All (Out & In)']);
    cb = cbAside(gca,'doy','k',[1 366]);
    set(gca,'colormap',colors);
    grid on;
    
    % frac out
    subplot(rows,cols,prc(cols,[iSex+1,3]));
    frac_out_cs = cumsum(frac_out_compiled')';
    seasonDoys = round(linspace(1,366,4+1));
    useDoys = [];
    for ii = 1:366
        if sum(frac_out_cs(ii,:)) > 0% && numel(unique(diff(odba_cs(ii,:)))) > 10 % remove some flat lines
            plot(frac_out_cs(ii,:),'color',[colors(ii,:) op],'linewidth',2);
            hold on;
        end
    end
    drawnow;
    xlim([1 size(frac_out_cs,2)]);
    ylim([0 15]);
    yticks([0,15]);
    xlabel('hours from sunrise');
    ylabel('cum sum frac out');
    set(gca,'fontsize',14);
    title([sexLabels{iSex+1},' ','Frac Out']);
    cb = cbAside(gca,'doy','k',[1 366]);
    set(gca,'colormap',colors);
    grid on;
end

% cum sum odba
subplot(rows,cols,prc(cols,[3,1]));
bothIds = sum(odba_compiled_M,2)~=0 & sum(odba_compiled_F,2)~=0;
x = find(bothIds);
lns = [];
ms = cumsum(sum(odba_compiled_M(bothIds,:),2));
ms_std = cumsum(mean(odba_compiled_M_std(bothIds,:),2));
lns(1) = plot(x,ms,'linewidth',2,'color',lineColors(1,:));
hold on;
%     plot(x,ms+ms_std,':','color',lineColors(1,:));
%     plot(x,ms-ms_std,':','color',lineColors(1,:));
fs = cumsum(sum(odba_compiled_F(bothIds,:),2));
fs_std = cumsum(mean(odba_compiled_M_std(bothIds,:),2));
lns(2) = plot(x,fs,'linewidth',2,'color',lineColors(2,:));
%     plot(x,fs+fs_std,':','color',lineColors(2,:));
%     plot(x,fs-fs_std,':','color',lineColors(2,:));
aa = plot(x,1,'k.','markersize',20);
lns(3) = aa(1);
ylabel('cum sum ODBA');

yyaxis right;
mm = mean(odba_compiled_M(bothIds,:),2);
plot(x,mm,'.','markersize',10,'color',lineColors(1,:));
hold on;
mf = mean(odba_compiled_F(bothIds,:),2);
plot(x,mf,'.','markersize',10,'color',lineColors(2,:));
ylabel('mean ODBA');
ylim([0 1]);
yticks(ylim);
set(gca,'ycolor','k')

xlabel('doy');
legend(lns(1:2),{'M','F'},'location','northwest');
legend box off;
title('annualized')
set(gca,'fontsize',14);
xlim([1 366]);

% cum sum asleep
subplot(rows,cols,prc(cols,[3,2]));
bothIds = sum(asleep_compiled_M,2)~=0 & sum(asleep_compiled_F,2)~=0;
x = find(bothIds);
lns = [];
ms = cumsum(sum(asleep_compiled_M(bothIds,:),2));
lns(1) = plot(x,ms,'linewidth',2,'color',lineColors(1,:));
hold on;
fs = cumsum(sum(asleep_compiled_F(bothIds,:),2));
lns(2) = plot(x,fs,'linewidth',2,'color',lineColors(2,:));
aa = plot(x,1,'k.','markersize',20);
lns(3) = aa(1);
ylabel('cum sum asleep');

yyaxis right;
mm = mean(asleep_compiled_M(bothIds,:),2);
plot(x,mm,'.','markersize',10,'color',lineColors(1,:));
hold on;
mf = mean(asleep_compiled_F(bothIds,:),2);
plot(x,mf,'.','markersize',10,'color',lineColors(2,:));
ylabel('mean asleep');
ylim([0 1]);
yticks(ylim);
set(gca,'ycolor','k')

xlabel('doy');
legend(lns(1:2),{'M','F'},'location','northwest');
legend box off;
title('annualized')
set(gca,'fontsize',14);
xlim([1 366]);


% fraction out of nest
subplot(rows,cols,prc(cols,[3,3]));
bothIds = sum(odba_compiled_M,2)~=0 & sum(odba_compiled_F,2)~=0;
x = find(bothIds);
lns = [];
ms = cumsum(mean(frac_out_compiled_M(bothIds,:),2));
lns(1) = plot(x,ms,'linewidth',2,'color',lineColors(1,:));
hold on;
fs = cumsum(mean(frac_out_compiled_F(bothIds,:),2));
lns(2) = plot(x,fs,'linewidth',2,'color',lineColors(2,:));
aa = plot(x,0,'k.','markersize',20);
lns(3) = aa(1);
ylabel('cum sum frac out');

yyaxis right;
mm = mean(frac_out_compiled_M(bothIds,:),2);
plot(x,mm,'.','markersize',10,'color',lineColors(1,:));
hold on;
mf = mean(frac_out_compiled_F(bothIds,:),2);
plot(x,mf,'.','markersize',10,'color',lineColors(2,:));
ylabel('mean frac out');
ylim([0 1]);
yticks(ylim);
set(gca,'ycolor','k')

xlabel('doy');
legend(lns(1:2),{'M','F'},'location','northwest');
legend box off;
title('annualized')
set(gca,'fontsize',14);
xlim([1 366]);