sqkey = readtable('sqkey.txt');
filespath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
load(fullfile(filespath,sqkey.filename{545}));
T = detect_sleepWake(T,2);
xs = 25552:27026;

op = 0.5;
colors = lines(5);
close all
ff(700,400);
plot(T.odba(xs),'k');
hold on;
plot(T.awake(xs)*.75-1.25,'color',[0 0 0 op]);
plot(strcmp(T.nest(xs),'Nest')*.75-2.5,'color',[1 0 0 op]);
ylim([-2.75 10]);
yyaxis right;
set(gca,'ycolor','r');
plot(T.temp(xs),'r');
ylim([9.5 30]);
xlim(size(xs));

xtickVals = linspace(1,numel(xs),12);
xtickLabelVals = circshift(0:2:23,-2);
xticks(xtickVals);
xticklabels(xtickLabelVals);
xlabel('hour of day');
yyaxis left;
yticks([-2.5,-1.75,-1.25,-.5,0:2:10]);
yticklabels({'Out Nest','In Nest','Asleep','Awake',0,2,4,6,8,10});
ylabel('ODBA (g)');

yyaxis right;
yticks(14:2:30);
ylabel('Collar Temp (C)');

set(gca,'fontsize',12);