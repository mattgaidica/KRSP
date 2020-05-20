nBins = 50;
close all
ff(600,600);
subplot(211);
histogram(meta.days_length/3600,nBins);
title(sprintf('day length (n=%i)',numel(meta.days_length)));
ylabel('count');
xlabel('time (h)');

subplot(212);
histogram(meta.nights_length/3600,nBins);
title(sprintf('night length (n=%i)',numel(meta.days_length)));
ylabel('count');
xlabel('time (h)');