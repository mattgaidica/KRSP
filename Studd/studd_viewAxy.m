% dataFile = '/Users/matt/Documents/Data/KRSP/doi_10.5061_dryad.1s1m8r7__v1/squirrelAxy_decisionTree.csv';
% A = readtable(dataFile);
close all
behaviors = unique(A.BEHAV);
contBehIds = unique(A.ID);

% split based on paper
nest = strcmp(A.BEHAV,'Nest');
notMoving = strcmp(A.BEHAV,'notMoving');
feeding = strcmp(A.BEHAV,'Feed');
foraging = strcmp(A.BEHAV,'CLIP') | strcmp(A.BEHAV,'DIGG') | strcmp(A.BEHAV,'CCONE') | strcmp(A.BEHAV,'SlowMove');
stationary = strcmp(A.BEHAV,'StatMove') | strcmp(A.BEHAV,'GRO') | strcmp(A.BEHAV,'vocal');
traveling = strcmp(A.BEHAV,'RunningMove');

classArr = {nest,notMoving,feeding,foraging,stationary,traveling};
classNames = {'nest','notMoving','feeding','foraging','stationary','traveling'};
ff(800,400);
colors = lines(numel(classArr));
for iC = 1:numel(classArr)
    scatter(A.tempC(classArr{iC}),A.ODBA(classArr{iC}),10,colors(iC,:),'filled');
    hold on;
end
legend(classNames);
% THIS IS ERROR PRONE, NEEDS CORRECTION FOR DAILY TEMPERATURE
figure;
anova1([A.tempC(feeding);A.tempC(nest)],[zeros(sum(feeding),1);ones(sum(nest),1)]);
anova1([A.tempC(foraging);A.tempC(nest)],[zeros(sum(foraging),1);ones(sum(nest),1)]);
anova1([A.tempC(traveling);A.tempC(nest)],[zeros(sum(traveling),1);ones(sum(nest),1)]);

ff(800,800);
scatter(A.tempC,A.ODBA,5,'filled');

ff(1200,300);
for iB = 1:numel(contBehIds)
    rows = A.ID == contBehIds(iB);
    ODBA = A.ODBA(rows);
    plot(ODBA);
    hold on;
end

% distribution of all behaviors
behCount = [];
for iB = 1:numel(behaviors)
    behCount(iB) = sum(strcmp(A.BEHAV,behaviors{iB}));
end
[v,k] = sort(behCount);
figure;
bar(behCount(k),'k');
xticks(1:numel(behaviors));
xticklabels(behaviors(k));
xtickangle(60);