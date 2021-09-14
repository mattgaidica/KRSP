function [T,W_z] = detect_sleepWake2(T,dls)
% dls is array of daylight seconds
% ex: dls = Tss.day_length(day(T.datetime,'dayofyear')) / 60; % min
if numel(dls) == 1
    error('use DLS version');
end
doPlot = true;
baseSec = 1440;
n = baseSec / 24;

dlsMinBlock = min([dls,baseSec-dls],[],2);
% is T.odba the mean for 60s or a decimated snapshot? !!it should do mean,
% more useful
% i.e. does it miss out on data, such that max is better?
% max could suffer from 'twitches' the skew the data
if doPlot
    colors = magma(n);
    close all;
    ff(900,900);
    subplot(211);
    plot(T.odba,'k');
    xlim([1 numel(T.odba)]);
%     xlim([1 3500]);
    ylim([0 12]);
    ylabel('OA');
    xlabel('Time (min)');
    set(gca,'fontsize',16);
    hold on;
    title(sprintf('%s - %s',datestr(T.datetime(1),'mmm dd YYYY'),datestr(T.datetime(end),'mmm dd YYYY')));
end
W = zeros(size(T.odba));
for iFilt = 1:n
%     disp(iFilt);
    filtFactor = 1440/iFilt;
    % this dynamically weighs the addition of the filter based on dls
    % *2 because dlsMinBlock is only half a 'cycle'
    thisSmooth = smoothdata(dlsMinBlock*2 >= filtFactor,'gaussian',baseSec/2)...
        .* smoothdata(T.odba,'loess',filtFactor);
    % must normalize
    W = W + normalize(thisSmooth,'range',[0 1]);
    if doPlot
        plot(W,'color',colors(iFilt,:));
    end
end
W_norm = normalize(W,'zscore'); % use this to estimate where sleep exists
useStd = std(W(W_norm < 0));
useMean = mean(W(W_norm < 0));
W_z = (W - useMean) ./ useStd;

if doPlot
    ylim([0 80]);
    yyaxis right;
    plot(xlim,[0,0],':k');
    hold on;
    plot(W_z,'k-','linewidth',1);
    plot(sign(W_z)+mean(ylim),'k-'); % binary sleep est
%     ylim([-1 3]);
    ylabel('Homeograph (H)');
    set(gca,'ycolor','k');
    
    subplot(212);
    histogram(T.odba(W_norm < 0),linspace(0,2,50));
    hold on;
    histogram(T.odba(W_norm >= 0),linspace(0,2,50));
    legend({'night OA','day OA'});
    xlabel('OA');
    ylabel('count');
    set(gca,'fontsize',16);
end

W_bin = zeros(numel(T.odba),1);
W_bin(W_z > 0) = 1;
T.odba_z = W_z;
T.awake = logical(W_bin);
T.asleep = ~W_bin;