load('algorithmApp_session.mat');
useRange = 1:numel(nest);

X = [];
X(:,1) = sense.tempZ(useRange);
X(:,2) = sense.tempGrad(useRange);
X(:,3) = sense.invOdba(useRange);
X(:,4) = sense.tempGradOdba(useRange);
X(:,5) = sense.tempRange(useRange);
kmeansNest = normalize(nest(useRange)','range',[0 1]);

% kmeans temp + odba
[idx,C] = kmeans(X(:,[1,3]),2);
tempOdbaNest = zeros(size(useRange));
tempOdbaNest(idx==1) = 1;

% pca
[coeff,score,latent,tsquared,explained] = pca(X);
disp(explained);
[idx,C] = kmeans(score(:,1:2),2); % use first two scores, usually captures most var
pcaNest = zeros(size(useRange));
pcaNest(idx==1) = 1;

% pca fixed
pcaFixedNest = fixTransitions(pcaNest',sense)';

close all;
ff(1200,400);
plotLabels = {'Orig','kmeansTempOdba','PCA','pcaFixedNest','binNestSense'};
plotData = {kmeansNest,tempOdbaNest,pcaNest,pcaFixedNest,binNestSense(useRange)'};
plotRange = 1:min([86400*3,numel(useRange)]);
yoff = 0:1.1:1.1*numel(plotData);
% fix kmeans dir
for ii = 2:numel(plotData) % compare to original nest
    if sum(plotData{1}==plotData{ii})/numel(plotData{1}) < 0.5
        plotData{ii} = ~plotData{ii};
    end
end
% plot and display corr
clc
for ii = 1:numel(plotData)
    plot(plotData{ii}(plotRange)+yoff(ii),'-','linewidth',2);
    hold on;
    for jj = 1:numel(plotData)
        str = sprintf("%s x %s = ",plotLabels{ii},plotLabels{jj});
        str = join(repmat(" ",[1,max(strlength(plotLabels))*2+7-strlength(str)]),"") + str;
        str = str + sprintf("%1.2f%%",100*sum(plotData{ii}==plotData{jj})/numel(plotData{ii}));
        fprintf("%s\n",str);
    end
    fprintf("\n");
end
op = 0.25;
transLocs = find(abs(diff(plotData{5}(plotRange)))==1);
for ii = 1:numel(transLocs)
    xline(transLocs(ii),'-','color',repmat(op,[1,4]));
end

legend(plotLabels,'autoupdate','off','location','northwest');
ylim([min(ylim)-1 max(ylim)+1]);
xlim([1 numel(plotRange)]);
yticklabels([]);

yyaxis right;
plot(odba(plotRange),'-','color',[0 0 1 op]);
hold on;
plot(normalize(X(plotRange,1),'range',[min(ylim),max(ylim)]),'-','color',[1 0 0 op],'linewidth',2);
set(gca,'ycolor','k');
set(gca,'fontsize',14);
yticklabels([]);


%%
[BF, BC] = bimodalitycoeff(alg_temp)
close all
ff(400,300);
histogram(alg_temp);

[BF, BC] = bimodalitycoeff(smOdba)
ff(400,300);
histogram(smOdba);

