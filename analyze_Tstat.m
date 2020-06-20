loadPath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
files = dir(fullfile(loadPath,'*.mat'));
nBins = (24*4); % 15min
binEdges = linspace(0,24,nBins+1) * 3600;
strFilts = {'Nest','Out'};
if do
    ees = zeros(2,nBins,366);
    ees_norm = zeros(2,366);
    ios = zeros(2,nBins,366);
    ios_norm = zeros(2,366);
    data_hist = cell(2,nBins,366);
    
    for iFile = 1:numel(files)
        disp(files(iFile).name);
        load(fullfile(loadPath,files(iFile).name));
        
        Tdatas = {T.odba,T.temp};
        for iData = 1:2
            doys = day(T.datetime,'dayofyear');
            for doy = unique(doys)'
                theseEntries = secDay(T.datetime(doys == doy));
                tdata = Tdatas{iData}(doys == doy);
                [~,~,bin] = histcounts(theseEntries,binEdges); % use counts for other than mean
                for iBin = 1:nBins
                    data_hist(iData,iBin,doy) = {[data_hist{iData,iBin,doy};tdata(bin == iBin)]};
                end
            end
        end
        
        for inOut = 1:2
            dts = Tstat.datetime(strcmp(Tstat.nest,strFilts{inOut})); % dates where Nest/Out exists
            doys = day(dts,'dayofyear');
            for doy = unique(doys)'
                ees_norm(inOut,doy) = ees_norm(inOut,doy) + 1;
                theseEntries = secDay(dts(doys == doy));
                counts = histcounts(theseEntries,binEdges);
                ees(inOut,:,doy) = ees(inOut,:,doy) + counts;
            end
            
            doys = day(Tstat.datetime,'dayofyear');
            for doy = unique(doys)'
                ios_norm(inOut,doy) = ios_norm(inOut,doy) + 1;
                doyDts = Tstat.datetime(doys == doy);
                theseEntries = secDay(doyDts);
                for iBin = 1:nBins
                    checkId = find(theseEntries <= binEdges(iBin+1),1,'last');
                    if isempty(checkId) % probably first few bins
                        if ~strcmp(strFilts{inOut},Tstat.nest{Tstat.datetime == doyDts(1)})
                            ios(inOut,iBin,doy) = ios(inOut,iBin,doy) + 1;
                        end
                    else
                        if strcmp(strFilts{inOut},Tstat.nest{Tstat.datetime == doyDts(checkId)})
                            ios(inOut,iBin,doy) = ios(inOut,iBin,doy) + 1;
                        end
                    end
                end
            end
        end
        
    end
    do = false;
    odba_hist = cellfun(@mean,squeeze(data_hist(1,:,:)));
    temp_hist = cellfun(@mean,squeeze(data_hist(2,:,:)));
end

Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2014.txt');

rows = 2;
cols = 3;
fs = 14;
titleStr = {'Enters','Exits';'In Nest','Out Nest'};
colors = lines(5);
op = 0.3;
close all;
ff(1600,800,2);
dataCols = {ees,ios};
normCols = {ees_norm,ios_norm};
Tcol = {temp_hist,odba_hist};
TcolTitle = {'Collar Temp (C)','Mean ODBA (g)'};

nTicks = 13;
for iCol = 1:2
    for inOut = 1:2
        data = squeeze(dataCols{iCol}(inOut,:,:)) ./ squeeze(normCols{iCol}(inOut,:));
        subplot(rows,cols,prc(cols,[inOut,iCol]));
        imagesc(data);
        set(gca,'ydir','normal');
        colormap(gca,parula);
        c = colorbar;
        c.Label.String = 'frac. of squirrels';
        caxis([0 1]);
        yticks(linspace(min(ylim),max(ylim),nTicks));
        hourOfDay = linspace(0,24,nTicks);
        yticklabels(datestr(datetime(1,1,1,floor(hourOfDay),...
            round((hourOfDay-floor(hourOfDay))*60),0),'HH:MM'));
        ylabel('hour of day');
        xlabel('doy');
        set(gca,'fontsize',fs);
        hold on;
        plot(secDay(Tss.sunrise)/(86400/nBins),'-','color',[1 1 1,op],'linewidth',2);
        plot(secDay(Tss.sunset)/(86400/nBins),'-','color',[1 1 1,op],'linewidth',2);
        title(titleStr{iCol,inOut});
    end
end
caxisVals = [20 30;0 1];
for iRow = 1:2
    subplot(rows,cols,prc(cols,[iRow,3]));
    data = Tcol{iRow};
    imagesc(data);
    set(gca,'ydir','normal');
    colormap(gca,parula);
    c = colorbar;
    c.Label.String = TcolTitle{iRow};
    caxis(caxisVals(iRow,:));
    ylabel('hour of day');
    xlabel('doy');
    set(gca,'fontsize',fs);
    hold on;
    plot(secDay(Tss.sunrise)/(86400/nBins),'-','color',[1 1 1,op],'linewidth',2);
    plot(secDay(Tss.sunset)/(86400/nBins),'-','color',[1 1 1,op],'linewidth',2);
    title(TcolTitle{iRow});
end