nDays = 7;
t = linspace(0,2*pi*nDays,nDays*24);

asleep_sim = square(t,70);
asleep_sim_naps = asleep_sim .* circshift(square(t,95),8);
close all;
ff(1200,600);

subplot(221);
plot(asleep_sim,'k','linewidth',3);
yticks(ylim);
yticklabels({'asleep','awake'});
xlabel('hours');
xticks(0:24:168);
xlim(size(asleep_sim));
grid on;
title('awake-asleep (no naps)');

subplot(223);
[c,lags] = xcorr(asleep_sim,'normalized');
plot(lags,c,'k','linewidth',3);
xlim([min(lags) max(lags)]);
xticks(-168:24:168);
xlabel('lag (hours)');
grid on
title('autocorrelation');

subplot(222);
plot(asleep_sim_naps,'k','linewidth',3);
yticks(ylim);
yticklabels({'asleep','awake'});
xlabel('hours');
xticks(0:24:168);
xlim(size(asleep_sim));
grid on;
title('awake-asleep (naps)');

subplot(224);
[c,lags] = xcorr(asleep_sim_naps,'normalized');
plot(lags,c,'k','linewidth',3);
xlim([min(lags) max(lags)]);
xticks(-168:24:168);
xlabel('lag (hours)');
grid on
title('autocorrelation');