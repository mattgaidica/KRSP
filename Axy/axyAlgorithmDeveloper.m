if do
    T = readtable('/Users/matt/Documents/Data/KRSP/AxyCSV/J1_Sep1_2014.csv');
    odba = T.odba;
    [temp,nest] = getTempAndNest(T.temp,60);
    zOffset = findZOffset(temp,odba,nest);
    do = 0;
end

wSweep = 0:0.5:5;
smSweep = [360,720,1200];

useThresh = 1;
useKmeans = 0;

if doCalc
    corrArr = [];
    lookupArr = [];
    ii = 0;
    for iW1 = 1:numel(wSweep)
        for iW2 = 1:numel(wSweep)
            for iW3 = 1:numel(wSweep)
                for iSm1 = 1:numel(smSweep)
                    for iSm2 = 1:numel(smSweep)
                        fprintf("%i%% -- %i/%i - %i/%i - %i/%i - %i/%i - %i/%i\n",...
                            round(100*ii/(numel(wSweep).^3*numel(smSweep).^2)),...
                            iW1,numel(wSweep),iW2,numel(wSweep),iW3,numel(wSweep),iSm1,numel(smSweep),iSm2,numel(smSweep));
                        wArr = [useKmeans,wSweep(iW1),wSweep(iW2),wSweep(iW3),smSweep(iSm1),smSweep(iSm2)];
                        [binNestSense,~,alg_temp,~,~,~,~,smOdba] = nestSenseAlg(temp,odba,nest,wArr,zOffset);
                        nMin = 20;
                        [t,periOdba,periTemp] = periEventNest(binNestSense,smOdba,alg_temp,nMin);
                        ii = ii + 1;
                        corrArr(ii) = mean(abs(mean(periOdba{1})-mean(periTemp{1}))) + mean(abs(mean(periTemp{2})-mean(periOdba{2})));
                        lookupArr(ii,1:5) = [iW1,iW2,iW3,iSm1,iSm2];
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