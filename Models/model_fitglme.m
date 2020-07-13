if do
    warning('off','all');
    sqkey = readtable('sqkey');
    filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
    Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
    Tss_doys = day(Tss.sunrise,'dayofyear');
    fitTable = table;
    iRow = 0;
    sqCount = 0;
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq})
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            T = detect_sleepWake(T,2);
            dtdoys = day(T.datetime,'dayofyear');
            undoys = unique(dtdoys);
            sqCount = sqCount + 1;
            for iDoy = 1:numel(undoys)
                dts = T.datetime(dtdoys==undoys(iDoy));
                [sunrise,sunset,day_length] = sunriseSunset(dts(1));
                % get weather, mean temperature, age
                % add previous day odba
                startId = closest(T.datetime,sunrise);
                if seconds(abs(T.datetime(startId)-sunrise)) < 5*60 % can not find sunrise
                    endId = startId + 1440 - 1; % use sunrise to sunrise as "day"
                    if endId <= numel(T.datetime)
                        iRow = iRow + 1;
                        nR = startId:endId;
                        fitTable.asleep(iRow) = sum(T.awake(nR))*60;
                        % this would be: does todays ODBA affect sleep?
                        % since the majority of night is included here
                        fitTable.odba(iRow) = sum(T.odba(nR));
                        maxVals = sort(T.odba_max(nR));
                        fitTable.odba_perf(iRow) = sum(maxVals(end-49:end));
                        % get next day odba: does todays sleep affect
                        % tomorrow's activity/performance?
                        if numel(T.odba) >= endId + 1440
                            fitTable.odba_next(iRow) = sum(T.odba(nR+1440));
                            maxVals = sort(T.odba_max(nR+1440));
                            fitTable.odba_next_perf(iRow) = sum(maxVals(end-49:end));
                        else
                            fitTable.odba_next(iRow) = NaN;
                            fitTable.odba_next_perf(iRow) = NaN;
                        end
                        
                        fitTable.doy(iRow) = undoys(iDoy);
                        fitTable.squirrelId(iRow) = {num2str(sqkey.squirrel_id(iSq))};
                        fitTable.dayLength(iRow) = day_length;
                        fitTable.sunrise(iRow) = secDay(sunrise);
                        fitTable.sunset(iRow) = secDay(sunset);
                        fitTable.sex(iRow) = sqkey.sex{iSq};
                    end
                end
            end
        end
    end
    warning('on','all');
    do = false;
end
clc
% for iRow = 1:size(fitTable,1)
%     fitTable.rnd(iRow) = normrnd(0,1);
% end

normType = 'zscore';
fitTable.asleep_norm = normalize(fitTable.asleep,normType);
fitTable.odba_norm = normalize(fitTable.odba,normType);
fitTable.odba_perf_norm = normalize(fitTable.odba_perf,normType);
fitTable.odba_next_norm = normalize(fitTable.odba_next,normType);
fitTable.odba_next_perf_norm = normalize(fitTable.odba_next_perf,normType);
fitTable.day_length_norm = normalize(fitTable.dayLength,normType);

close all
ff(1200,1200);
rows = 2;
cols = 3;
subplot(rows,cols,1);
histogram(fitTable.asleep_norm);
title('asleep_norm');
subplot(rows,cols,2);
histogram(fitTable.odba_norm);
title('odba_norm');
subplot(rows,cols,3);
histogram(fitTable.odba_perf_norm);
title('odba_perf_norm');
subplot(rows,cols,4);
histogram(fitTable.odba_next_norm);
title('odba_next_norm');
subplot(rows,cols,5);
histogram(fitTable.odba_next_perf_norm);
title('odba_next_perf_norm');
subplot(rows,cols,6);
histogram(fitTable.day_length_norm);
title('day_length_norm');

glme = fitglme(fitTable,'asleep_norm ~ odba_norm^2 + odba_next_norm^2 + odba_perf_norm + odba_next_perf_norm + (1|sex)',...
    'Distribution','Normal','DummyVarCoding','effects')
[psi,dispersion,stats] = covarianceParameters(glme);

% close all
ff(1200,400);
subplot(131);
plotResiduals(glme,'histogram');
subplot(132);
plotResiduals(glme,'fitted');
subplot(133);
plotResiduals(glme,'probability');

stats{1}
    
%     doy + squirrelId + dayLength + sunrise + sunset + (1|sex)