if do
    sqkey = readtable('sqkey');
    filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
    Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
    Tss_doys = day(Tss.sunrise,'dayofyear');
    
    nBins = 24; % 15min
    binEdges = linspace(0,24,nBins+1) * 3600;
    
    doy_arr = [];
    odba_arr = zeros(1,nBins);
    odba_arr_out = zeros(1,nBins);
    odba_arr_nest = zeros(1,nBins);
    sqCount = 0;
    doyCount = 0;
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq})
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            if isValidT(T,true)
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
                    theseSecs = secDay(dts);
                    thisSunrise = secDay(Tss.sunrise(Tss_doys==undoys(iDoy)));
                    [v,k] = sort(abs(theseSecs-thisSunrise));
                    sunriseOdba = theseOdba(k(1):k(2));
                    sunriseSecs = theseSecs(k(1):k(2));
                    sunriseSecs = sunriseSecs - sunriseSecs(1);
                    addId = find(sunriseSecs<0,1,'first');
                    sunriseSecs(addId:end) = sunriseSecs(addId:end) + sunriseSecs(addId-1) - sunriseSecs(addId);
                    sunriseNest = theseNest(k(1):k(2));
                    for iBin = 1:nBins
                        odba_arr(doyCount,iBin) = mean(sunriseOdba(...
                            sunriseSecs >= binEdges(iBin) & sunriseSecs < binEdges(iBin+1)));
                        odba_arr_out(doyCount,iBin) = mean(sunriseOdba(...
                            sunriseSecs >= binEdges(iBin) & sunriseSecs < binEdges(iBin+1) &...
                            strcmp(sunriseNest,'Out')));
                        odba_arr_nest(doyCount,iBin) = mean(sunriseOdba(...
                            sunriseSecs >= binEdges(iBin) & sunriseSecs < binEdges(iBin+1) &...
                            strcmp(sunriseNest,'Nest')));
                    end
                    doy_arr(doyCount) = undoys(iDoy);
                end
                
                % %                 for doy = unique(dtdoys)'
                % %                     doyCount = doyCount + 1;
                % %                     theseDts = secDay(T.datetime(dtdoys==doy));
                % %                     theseOdba = T.odba(dtdoys==doy);
                % %                     theseNest = T.nest(dtdoys==doy);
                % %                     for iBin = 1:nBins
                % %                         odba_arr(doyCount,iBin) = mean(theseOdba(...
                % %                             theseDts >= binEdges(iBin) & theseDts < binEdges(iBin+1)));
                % %                         odba_arr_out(doyCount,iBin) = mean(theseOdba(...
                % %                             theseDts >= binEdges(iBin) & theseDts < binEdges(iBin+1) &...
                % %                             strcmp(theseNest,'Out')));
                % %                         odba_arr_nest(doyCount,iBin) = mean(theseOdba(...
                % %                             theseDts >= binEdges(iBin) & theseDts < binEdges(iBin+1) &...
                % %                             strcmp(theseNest,'Nest')));
                % %                     end
                % %                     doy_arr(doyCount) = doy;
                % %                 end
            end
        end
    end
    do = false;
end

close all
op = 0.1;
data_labels = {'All','Out of Nest','In Nest'};
data_sets = {odba_arr,odba_arr_out,odba_arr_nest};
colors = [parula(ceil(366/2));flip(parula(floor(366/2)))];
ylimVals = [12,18,4];
ff(1200,500);
for iData = 1:3
    odba_compiled = zeros(366,nBins);
    for ii = 1:366
        odba_compiled(ii,:) = nanmean(data_sets{iData}(doy_arr==ii,:));
    end
    
    subplot(1,3,iData);
    odba_cs = cumsum(odba_compiled')';
    for ii = 1:366
        if sum(odba_cs(ii,:)) > 0
            plot(odba_cs(ii,:),'color',[colors(ii,:) op],'linewidth',2);
            hold on;
            drawnow;
        end
    end
    xlim([1 size(odba_cs,2)]);
    ylim([0 ylimVals(iData)]);
    yticks([0:2:ylimVals(iData)]);
    xlabel('hours from sunrise');
    ylabel('cum. sum. ODBA');
    set(gca,'fontsize',14);
    title(data_labels{iData});
    grid on;
end