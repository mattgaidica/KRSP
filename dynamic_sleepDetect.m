%% test data

T = table;
T.odba = zeros(1440*3,1);
T.odba(361:1080) = 1;
T.odba((361:1080)+1440) = 1;
T.odba((301:500)+1440*2) = 1;

dls = [repmat(numel(361:1080),[1440,1]);repmat(numel(361:1080),[1440,1]);repmat(numel(301:500),[1440,1])];
dlsMinBlock = min([dls,1440-dls],[],2);

close all
ff(1400,900);
plot(T.odba,'k-','linewidth',2);
a = normalize(smoothdata(T.odba,'gaussian',1440));
% a = a - movmean(a,1440);
yyaxis right;
plot(a);

%% use on real data
dls = Tss.day_length(day(T.datetime,'dayofyear')) / 60; % min
dlsMinBlock = min([dls,1440-dls],[],2);

n = 60;

colors = jet(n);
close all;
ff(1400,900);
axs = [];
axs(1) = subplot(211);
plot(T.odba,'k');
% xlim([1 3500]);
xlim([1 numel(T.odba)]);
ylim([0 12]);
ylabel('\DeltaOA');
xlabel('Time (min)');
set(gca,'fontsize',14);
hold on;

W = zeros(size(T.odba));
for iFilt = 1:n
    filtFactor = 1440/iFilt;
    thisSmooth = smoothdata(dlsMinBlock*2 >= filtFactor,'gaussian',720) .* smoothdata(T.odba,'gaussian',filtFactor);
    W = W + normalize(thisSmooth,'range',[0,1]);
    plot(W,'color',colors(iFilt,:));
end
ylim auto;
W_norm = normalize(W,'zscore'); % use this to estimate where sleep exists
useStd = std(W(W_norm < 0));
useMean = mean(W(W_norm < 0));
W_z = (W - useMean) ./ useStd;

%     ylim([0 80]);
yyaxis right;
plot(xlim,[0,0],':k');
hold on;
plot(W_z,'k-','linewidth',2);
plot(sign(W_z)+mean(ylim),'k-'); % binary sleep est
%     ylim([-1 3]);
ylabel('homeograph Z-score');
set(gca,'ycolor','k');

axs(2) = subplot(212);
plot(T.odba,'k');
xlim([1 numel(T.odba)]);
linkaxes(axs,'x');