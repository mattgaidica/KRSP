function [T,W_z] = detect_sleepWake2(T,n)
doPlot = false;

% is T.odba the mean for 60s or a decimated snapshot? !!it should do mean,
% more useful
% i.e. does it miss out on data, such that max is better?
% max could suffer from 'twitches' the skew the data
colors = jet(n);
if doPlot
    close all;
    ff(1200,900);
    subplot(211);
    plot(T.odba,'k');
    xlim([1 3500]);
    ylim([0 12]);
    ylabel('\DeltaOA');
    xlabel('Time (min)');
    set(gca,'fontsize',14);
    hold on;
end
W = zeros(size(T.odba));
for ii = 1:n
    W = W + smoothdata(T.odba,'loess',1440/ii);
    if doPlot
        plot(W,'color',colors(ii,:));
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
%     ylim([-1 3]);
    ylabel('homeograph Z-score');
    set(gca,'ycolor','k');
    subplot(212);
    histogram(T.odba(W_norm < 0),100);
    hold on;
    histogram(T.odba(W_norm >= 0),100);
    legend({'night ODBAs','day ODBAs'});
end

W_bin = zeros(numel(T.odba),1);
W_bin(W_z > 0) = 1;
T.odba_z = W_z;
T.awake = logical(W_bin);
T.asleep = ~W_bin;