%% setup predictor data
if do
    T = readtable('/Volumes/GAIDICASSD/KRSP/KRSP Axy Data/Emily Studds Axy/April 2015 - Week 1/B10_Apr23_2015.csv');
    load('trainingNestApril2015.mat');
    do = 0;
    odba = T.odba;
    [temp,nest] = getTempAndNest(T.temp,60);
end
nInt = 10;
useRange = 1:nInt:numel(temp);
Y = trainingNest(useRange);
X = []; % numel(nest) rows, cols are predictors
ii = 1;

predictorLabels = {};

X(:,ii) = normalize(temp(useRange));
predictorLabels{ii} = "temp";
ii = ii+1;

smoothVals = [30,120,600]/nInt;% 1:120:60*10;
for jj = smoothVals
    smGrad = gradient(smoothdata(temp(useRange),'gaussian',jj));
    X(:,ii) = normalize(smGrad);
    predictorLabels{ii} = sprintf("tempGrad-sm%i",jj*nInt);
    ii = ii+1;
end

for jj = smoothVals
    smGrad = gradient(smoothdata(temp(useRange),'gaussian',jj));
    X(:,ii) = normalize(smGrad.*abs(smGrad));
    predictorLabels{ii} = sprintf("tempGrad^2-sm%i",jj*nInt);
    ii = ii+1;
end

smoothVals = [120,600,1200,2400]/nInt;% 1:120:60*10;
for jj = smoothVals
    X(:,ii) = -smoothdata(normalize(odba(useRange)),'gaussian',jj);
    predictorLabels{ii} = sprintf("ODBA-sm%i",jj*nInt);
    ii = ii+1;
end
%% fit model
SVMModel = fitcsvm(X,Y);
beep;
%% ^plot betas
% close all;
ff(800,400);
betas = SVMModel.Beta;
bar(1:numel(betas),betas,'facecolor','k');
xticks(1:numel(betas));
xticklabels(predictorLabels);
xtickangle(30);
ylabel('Beta');
set(gca,'fontsize',14);