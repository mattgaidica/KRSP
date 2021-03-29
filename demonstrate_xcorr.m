nDays = 7;
t = linspace(0,2*pi*nDays,nDays*1440);

asleep_sim = square(t,60);
asleep_sim_naps = asleep_sim .* -circshift(square(t,5),400);
% close all;
ff(1200,600);

subplot(221);
plot(asleep_sim,'k','linewidth',3);
yticks(ylim);
yticklabels({'asleep','awake'});
xlabel('hours');
xticks(0:1440:1440*nDays);
xlim(size(asleep_sim));
grid on;
title('awake-asleep (no naps)');

subplot(223);
[c,lags] = xcorr(asleep_sim,'normalized');
plot(lags,c,'k','linewidth',3);
xlim([0 1440*3]);
xticks(-1440*nDays:1440:1440*nDays);
xlabel('lag (hours)');
grid on
title('autocorrelation');

subplot(222);
plot(asleep_sim_naps,'k','linewidth',3);
yticks(ylim);
yticklabels({'asleep','awake'});
xlabel('hours');
xticks(0:1440:1440*nDays);
xlim(size(asleep_sim));
grid on;
title('awake-asleep (naps)');

subplot(224);
[c,lags] = xcorr(asleep_sim_naps,'normalized');
plot(lags,c,'k','linewidth',3);
xlim([0 1440*3]);
xticks(-1440*nDays:1440:1440*nDays);
xlabel('lag (hours)');
grid on
title('autocorrelation');