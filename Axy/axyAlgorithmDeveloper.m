if do
    T = readtable('/Users/matt/Documents/Data/KRSP/AxyCSV/J1_Sep1_2014.csv');
    odba = T.odba;
    [temp,nest] = getTempAndNest(T.temp,60);
    offset_odba = findZOffset(temp,odba,nest);
    do = 0;
end

threshSweep = .2;
wSweep = 0:1:5;
smSweep = [360,720,1200];

if doCalc
    corrArr = [];
    lookupArr = [];
    ii = 0;
    for iT = 1:numel(threshSweep)
        for iW1 = 1:numel(wSweep)
            for iW2 = 1:numel(wSweep)
                for iW3 = 1:numel(wSweep)
                    for iSm1 = 1:numel(smSweep)
                        for iSm2 = 1:numel(smSweep)
                            fprintf("%i%% -- %i/%i - %i/%i - %i/%i - %i/%i - %i/%i - %i/%i\n",...
                                round(100*ii/(numel(threshSweep)*numel(wSweep).^3*numel(smSweep).^2)),...
                                iT,numel(threshSweep),iW1,numel(wSweep),iW2,numel(wSweep),iW3,numel(wSweep),iSm1,numel(smSweep),iSm2,numel(smSweep));
                            
                            senseParams = {};
                            senseParams.thresh = threshSweep(iT);
                            senseParams.w_kmeans = 0;
                            senseParams.w_temp = wSweep(iW1);
                            senseParams.w_tempGrad = wSweep(iW2);
                            senseParams.w_odba = wSweep(iW3);
                            senseParams.sm_tempGrad = smSweep(iSm1);
                            senseParams.sm_odba = smSweep(iSm2);
                            senseParams.offset_odba = offset_odba;
                            [binNestSense,sense] = nestSenseAlg(temp,odba,nest,senseParams);
                            
                            ii = ii + 1;
                            corrArr(ii) = sum(binNestSense==nestRange) / numel(nestRange);
                            lookupArr(ii,1:6) = [iT,iW1,iW2,iW3,iSm1,iSm2];
                        end
                    end
                end
            end
        end
    end
    doCalc = 0;
end

%%
[y,k] = sort(corrArr);
maxk = k(end);
fprintf("Best Corr (r=%1.4f): w_Temp:%1.1f, w_TempGrad:%1.1f, w_ODBA:%1.1f, smGrad:%i, smODBA:%i\n",...
    y(end),wSweep(lookupArr(maxk,1)),wSweep(lookupArr(maxk,2)),wSweep(lookupArr(maxk,3)),...
    smSweep(lookupArr(maxk,4)),smSweep(lookupArr(maxk,5)));

close all
ff(1200,300);
titleLabels = {'Temp (w)','Temp Gradient (w)','ODBA (w)','smGrad','smODBA'};
ylabelLabel = 'Corr w/ Studd';
ms = 7;
for ii = 1:5
    subplot(1,5,ii);
    if ii < 4
        x = wSweep(lookupArr(k,ii));
    else
        x = smSweep(lookupArr(k,ii));
    end
    scatter(x,y,ms,magma(numel(y)),'filled');
    xticks(unique(x));
    xlabel(titleLabels{ii});
    ylabel(ylabelLabel);
    title(titleLabels{ii});
    grid on;
    set(gca,'fontsize',14);
    %     ylim([0 1]);
end
% saveas(gcf,fullfile(savePath,'model-corr-analysis.jpg'));