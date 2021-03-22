% setup with /Users/matt/Documents/MATLAB/KRSP/predict_awake.m
close all
% statsMult = [1 -1;-1 1;1 1;-1 -1];
statsMult = [0 -1;0 1;1 0;-1 0];
ff(800,400);

subplot(121);
op = 0.075;
useIds = [2,3,1,4,2];
colors = mycmap('/Users/matt/Documents/MATLAB/KRSP/util/seasons2.png',5);

lns = [];
for ii = 1:size(overlapStats,1)
    theseStats = overlapStats(ii,:);
    for jj = 1:4
        xs = [theseStats(useIds(jj))*statsMult(useIds(jj),1),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),1)];
        ys = [theseStats(useIds(jj))*statsMult(useIds(jj),2),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),2)];
        lns(1) = plot(xs,ys,'color',[0,0,0,op]);
        hold on;
    end
end

sIds = round(linspace(1,366,5));
seasonDoys = circshift(1:366,60);
for iSeason = 1:4
    ss = ismember(mean_doys,seasonDoys(sIds(iSeason):sIds(iSeason+1)));
    theseStats = mean(overlapStats(ss,:));
    for jj = 1:4
        xs = [theseStats(useIds(jj))*statsMult(useIds(jj),1),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),1)];
        ys = [theseStats(useIds(jj))*statsMult(useIds(jj),2),theseStats(useIds(jj+1))*statsMult(useIds(jj+1),2)];
        lns(iSeason+1) = plot(xs,ys,'color',[colors(iSeason,:),1],'linewidth',2);
        hold on;
    end
end
legend(lns,{'Individual','Winter','Spring','Summer','Autumn'});
legend box off
title('Nest-sleep Overlap');

xlim([-0.5 0.75]);
ylim([-0.75 0.75]);
xticks(sort([0,xlim]));
yticks(xticks);
xticklabels(abs(xticks));
yticklabels(xticklabels);
offset = 0.05;
fs = 14;
text(0,max(ylim)-offset*2,'in-asleep','horizontalalignment','center','fontsize',fs);
text(max(xlim),-offset*2,'out-awake','horizontalalignment','right','fontsize',fs);
text(0,min(ylim)+offset*2,'in-awake','horizontalalignment','center','fontsize',fs);
text(min(xlim)+offset,offset*2,'out-asleep','horizontalalignment','left','fontsize',fs);
% grid on;
xlabel('fraction of day');
ylabel('fraction of day');
set(gca,'fontsize',16);
box off;

subplot(122);
for iSeason = 1:4
    ss = ismember(mean_doys,seasonDoys(sIds(iSeason):sIds(iSeason+1)));
    x = sq_inNestMin(ss)*24;
    y = sq_asleepMin(ss)*24;
    plot(x,y,'.','color',colors(iSeason,:),'markersize',15);
    hold on;
end
xlim([0 24]);
xticks(0:4:24);
ylim([2 16]);
f = fit(sq_inNestMin'*24,sq_asleepMin'*24,'poly1');
% f(x) = p1*x + p2
lns = [];
lns(1) = plot(NaN,NaN,'k.','markersize',15);
lns(2) = plot(xlim,[f.p1*min(xlim)+f.p2,f.p1*max(xlim)+f.p2],'k-');
ci = confint(f);
lns(3) = plot(xlim,[ci(1,1)*min(xlim)+ci(1,2),ci(1,1)*max(xlim)+ci(1,2)],'k:');
plot(xlim,[ci(2,1)*min(xlim)+ci(2,2),ci(2,1)*max(xlim)+ci(2,2)],'k:');

[r,p] = corr(sq_inNestMin'*24,sq_asleepMin'*24);
title('Asleep vs. In Nest');
xlabel('In Nest (hrs/day)');
ylabel('Asleep (hrs/day)');
grid on;
set(gca,'fontsize',16);
legend(lns,{'Squirrel Mean','Linear Fit','95% Confidence'},'location','northwest');
text(12,3,sprintf('r = %1.2f, p = %1.2e',r,p),'horizontalalignment','center');
