filespath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
if do
    files = dir(fullfile(filespath,'*.mat'));
    overlapStats = [];
    iCount = 0;
    mean_doys = [];
    sq_inNestMin = [];
    sq_asleepMin = [];
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq})
            load(fullfile(loadPath,sqkey.filename{iSq}));
            if ~sqkey.isValid(iSq) || ~isValidT(T,true) % check temp
                disp('invalid');
                continue;
            end
        else
            disp('not loading');
            continue;
        end
        fprintf("%i/%i - %s\n",iSq,size(sqkey,1),sqkey.filename{iSq});
        iCount = iCount + 1;
        T.datetime = T.datetime + minutes(sqkey.shiftMin(iSq));
        dls = Tss.day_length(day(T.datetime,'dayofyear')) / 60; % min
        T = detect_sleepWake2(T,dls);
        
        T.nest_bin = strcmp(T.nest,'Nest');
        overlapStats(iCount,1) = sum(T.nest_bin & T.awake) / size(T,1); % in-awake
        overlapStats(iCount,2) = sum(T.nest_bin & ~T.awake) / size(T,1); % in-asleep
        overlapStats(iCount,3) = sum(~T.nest_bin & T.awake) / size(T,1); % out-awake
        overlapStats(iCount,4) = sum(~T.nest_bin & ~T.awake) / size(T,1); % out-asleep
        mean_doys(iCount) = mean(unique(day(T.datetime,'dayofyear')));
        sq_inNestMin(iCount) = sum(T.nest_bin) ./ size(T,1);
        sq_asleepMin(iCount) = sum(T.asleep) ./ size(T,1);
    end
    do = false;
end

if false
    bee_stats = overlapStats(:);
    bee_cats = repmat(1:4,[size(bee_stats,1)/4,1]);
    bee_cats = bee_cats(:);
    catLabels = {'In & Awake','In & Asleep','Out & Awake','Out & Asleep'};
    
    close all
    ff(500,500);
    beeswarm(bee_cats,bee_stats,'corral_style','omit','overlay_style','sd');
    axis tight
    xticks([1:4]);
    xticklabels(catLabels);
    ylabel('frac. of time');
end

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
title(sprintf('r = %1.2f, p = %1.2e',r,p));
xlabel('In Nest (hrs/day)');
ylabel('Asleep (hrs/day)');
grid on;
set(gca,'fontsize',16);
legend(lns,{'Squirrel Mean','Linear Fit','95% Confidence'},'location','northwest');
