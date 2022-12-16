% function findTempDelay(T)

searchWindow = 60*15; % minutes
nest = strcmp(T.Nest,'Nest');
odba = T.odba;
% temp = T.temp;
% temp_filt = temp;
% temp_filt(temp==22.5) = NaN;
% temp_filt = inpaint_nans(temp_filt);
% temp_smooth = smoothdata(temp_filt,'gaussian',60);

sumODBA = [];
for ii = 1:searchWindow
    useRange = 1:numel(nest)-searchWindow;
    nestShift = circshift(nest,-ii);
    nestSelect = nestShift(useRange);
    odbaSelect = odba(useRange);
    sumODBA(ii) = mean(odbaSelect(nestSelect==1));
end
%% ^
t = linspace(0,searchWindow/60,searchWindow);
close all
ff(800,300);
plot(t,sumODBA,'k-','linewidth',2);
hold on;
[v,k] = min(sumODBA);
plot(t(k),v,'x','markerSize',20);
set(gca,'fontsize',14);
xlabel('Lag (min)');
grid on;
title(sprintf('Minimizing ODBA In Nest - ODBA lags Nest by %1.0f sec / %1.2f min',k,k/60));
ylabel('ODBA');
% saveas(gcf,fullfile(savePath,'nest-axy-lag.jpg'));