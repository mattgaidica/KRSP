% !! it seems that you want to check the data around a transition, not just
% at the transition
threshSweep = .5;
wSweep = 0:2:4;
smSweep = 60:360:360*4;

nestRange = trainingNest;
odbaRange = odba;
tempRange = temp;
% odbaRange = odba;
% tempRange = temp;
if do
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
                            binNestSense = nestSenseAlg(tempRange,odbaRange,...
                                [threshSweep(iT),wSweep(iW1),wSweep(iW2),wSweep(iW3),smSweep(iSm1),smSweep(iSm2)]);
                            ii = ii + 1;
                            corrArr(ii) = sum(binNestSense==nestRange) / numel(nestRange);
                            lookupArr(ii,1:6) = [iT,iW1,iW2,iW3,iSm1,iSm2];
                        end
                    end
                end
            end
        end
    end
    do = false;
end
%Corr (r=0.9812): Thresh=0.5, w_Temp:8.0, w_TempGrad:4.0, w_ODBA:8.0, smGrad:560, smODBA:360
%%
[y,k] = sort(corrArr);
maxk = k(end);
fprintf("Best Corr (r=%1.4f): Thresh=%1.1f, w_Temp:%1.1f, w_TempGrad:%1.1f, w_ODBA:%1.1f, smGrad:%i, smODBA:%i\n",...
    y(end),threshSweep(lookupArr(maxk,1)),wSweep(lookupArr(maxk,2)),wSweep(lookupArr(maxk,3)),wSweep(lookupArr(maxk,4)),...
    smSweep(lookupArr(maxk,5)),smSweep(lookupArr(maxk,6)));

close all
ff(1200,300);
titleLabels = {'Threshold (Z)','Temp (w)','Temp Gradient (w)','ODBA (w)','smGrad','smODBA'};
ylabelLabel = 'Corr w/ Studd';
ms = 7;
for ii = 1:6
    subplot(1,6,ii);
    if ii == 1
        x = threshSweep(lookupArr(k,ii));
    elseif ii < 5
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
    ylim([.9 1]);
end
% saveas(gcf,fullfile(savePath,'model-corr-analysis.jpg'));
%%
startSample = 1+60*60* 48;
useSamples = 60*60*48; % hours
useRange = startSample:startSample+useSamples-1;
tempRange = temp(useRange);
oldNestRange = strcmp(T.Nest(useRange),'Nest');
odbaRange = odba(useRange);
% nestRange = nest(useRange);
nestRange = strcmp(T.Nest(useRange),'Nest'); % compare against original
% % ax = [];
% % colors = lines(5);
% % close all
% % ff(1200,800);
% % ax(1) = subplot(211);
% % plot(odbaRange,'k');
% % hold on;
% % plot(normalize(gradient(tempRange)),'r-');
% % plot(normalize(tempRange),'color',[1 0 0 0.5],'linewidth',2);
% % yyaxis right;

[binNestSense,nestSense] = nestSenseAlg(tempRange,odbaRange,[]);

% % binNestSense = normalize(sign(nestSense),'range');
% % binNestSense(binNestSense==0.5) = NaN;
% % binNestSense = fillmissing(binNestSense,'previous');
ax = [];
t = linspace(0,numel(nestSense)/60,numel(nestSense));
close all
ff(1200,900);
ax(1) = subplot(211);
ln1 = plot(t,odbaRange,'-');
hold on;
ln2 = plot(t,nestSense,'k-','linewidth',3);
ln3 = yline(useThresh,'k:','linewidth',1);
grid on;
set(gca,'fontsize',14);
xlabel('Time (min)');
ylabel('p(Nest)');

yyaxis right;
ln4 = plot(t,normalize(tempRange),'-','color','r','linewidth',1);
hold on;
ln5 = plot(t,binNestSense,'-','color',colors(5,:),'linewidth',2);
ln6 = plot(t,nestRange,'-','color',[colors(4,:),0.4],'linewidth',2);
set(gca,'ycolor','k');
yticks([0,1]);
ylim([-4 2]);
yticklabels({'Out of Nest','In Nest'});
legend([ln1,ln2,ln4,ln5,ln6],{'ODBA','p(Nest)','Temp','New Nest','Old Nest'});

ax(2) = subplot(212);
plot(t,comp1,'-');
hold on;
plot(t,comp2,'-');
plot(t,comp3,'-');
plot(t,nestSense,'k-','linewidth',3);
grid on;
legend({'Temp','Grad','ODBA','Sense'});
set(gca,'fontsize',14);
yline(0,'k:');

linkaxes(ax,'x');
% xlim(usexlim);
% saveas(gcf,fullfile(savePath,'algorithmic-threshold-example.jpg'));