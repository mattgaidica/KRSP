%%
close all
ff(1500,900);

colors = lines(5);
% nR = 1:1440*7;
% data = smoothdata(T.odba_max,'gaussian',60);
data = T.odba;
n = 60;

subplot(211);
plot(normalize(data),'color',repmat(0.2,[1,4]));
hold on;
tryMethods = {'loess'};

for jj = 1:numel(tryMethods)
    d1 = zeros(size(data));
    for ii = 1:4:n
        d1 = d1 + smoothdata(data,tryMethods{jj},1440/ii);
    end
    zd = (d1-mean(d1))./std(d1);
    plot(zd,'linewidth',1.5);
    drawnow;
end
legend([{'axy'},tryMethods(:)']);

subplot(212);
plot(zd,'color',colors(1,:),'linewidth',1.5);
hold on;
xlim([1 size(zd,1)]);
plot([1 size(zd,1)],[0 0],'k:','linewidth',1);
ylabel('homeogram Z-score');

yyaxis right;

dayAct = data;
dayAct(zd < 0) = NaN;

nightAct = data;
nightAct(zd >= 0) = NaN;

plot(dayAct,'k-');
hold on;
plot(nightAct,'r-');
grid on;
ylabel('axy amplitude');
set(gca,'ycolor','k');

xlabel('Time (min)');

legend({'homeograph','<0 ~ asleep','axy "awake"','axy "asleep"'},'location','northwest');
set(gca,'fontsize',16);

% figure;
% histogram(dayAct,linspace(0,3,100));
% hold on;
% histogram(nightAct,linspace(0,3,100));
% set(gca,'yscale','log');


%%

ff(1200,600);

subplot(211);
plot(data,'k','linewidth',1.5);
yyaxis right;
gsmooth = smoothdata(data,'gaussian',1440);
gz = (gsmooth - mean(gsmooth)) ./ std(gsmooth);
mvme = movmean(gz,1440);
plot(gz,'color',colors(1,:),'linewidth',1.5);
hold on;
plot(mvme,'--','color',colors(2,:),'linewidth',1.5);
xlim(size(nR));
grid on;

adjGz = gz - mvme;


subplot(212);
plot(adjGz,'color',colors(1,:),'linewidth',1.5);
yyaxis right;

dayAct = data;
dayAct(adjGz < 0) = NaN;

nightAct = data;
nightAct(adjGz >= 0) = NaN;

plot(dayAct,'k');
hold on;
plot(-nightAct,'r');
% xlim(size(nR));
grid on;

figure;
histogram(dayAct,1:0.1:10);
hold on;
histogram(nightAct,1:0.1:10);
set(gca,'yscale','log');