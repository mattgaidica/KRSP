if do
    sqkey = readtable('sqkey');
    filePath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
    % constant datetime x-axis?
    
    xs = zeros(size(sqkey,1),1);
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq})
            disp(sqkey.filename{iSq});
            load(fullfile(filePath,sqkey.filename{iSq})); % T, Tstat
            xs(iSq) = day(T.datetime(1),'dayofyear');
        end
    end
    [v,k] = sort(xs);

    M_scatt_x = [];
    M_scatt_y = [];
    F_scatt_x = [];
    F_scatt_y = [];
    F_scattLac_x = [];
    F_scattLac_y = [];
    F_scattPreg_x = [];
    F_scattPreg_y = [];
    sqCount = 0;
    for iSq = 1:size(sqkey,1)
        useRow = k(iSq);
        if ~isempty(sqkey.filename{useRow})
            disp(sqkey.filename{iSq});
            sqCount = sqCount + 1;
            load(fullfile(filePath,sqkey.filename{useRow})); % T, Tstat
            doys = unique(day(T.datetime,'dayofyear'));

            if strcmp(sqkey.sex{useRow},'M')
                M_scatt_x = [M_scatt_x;doys];
                M_scatt_y = [M_scatt_y;repmat(sqCount,[numel(doys),1])];
            else
                F_scatt_x = [F_scatt_x;doys];
                F_scatt_y = [F_scatt_y;repmat(sqCount,[numel(doys),1])];
                if strcmp(sqkey.sex_status{useRow},'lactating')
                    F_scattLac_x = [F_scattLac_x;doys];
                    F_scattLac_y = [F_scattLac_y;repmat(sqCount,[numel(doys),1])];
                elseif strcmp(sqkey.sex_status{useRow},'pregnant')
                    F_scattPreg_x = [F_scattPreg_x;doys];
                    F_scattPreg_y = [F_scattPreg_y;repmat(sqCount,[numel(doys),1])];
                end
            end
        end
    end
    do = false;
end


close all;
ff(1400,900,2);
subplot(211);
ms = 5;
plot(F_scatt_x,F_scatt_y,'m.','markersize',ms);
hold on;
plot(M_scatt_x,M_scatt_y,'b.','markersize',ms);
ylim([1 sqCount]);
xlim([1 366]);
xlabel('doy');
ylabel('squirrel');
set(gca,'fontsize',14);
title(sprintf('%i squirrels (%iF + %iM), %i rec days',sqCount,...
    numel(unique(F_scatt_y)),numel(unique(M_scatt_y)),numel(M_scatt_x)+numel(F_scatt_y)));

lw = 2;
subplot(212);
binEdges = linspace(0.5,366.5,366+1);
nM = histcounts(M_scatt_x,binEdges);
nF = histcounts(F_scatt_x,binEdges);
nFLac = histcounts(F_scattLac_x,binEdges);
nFPreg = histcounts(F_scattPreg_x,binEdges);
plot(1:366,nM+nF,'-k','linewidth',lw);
hold on;
plot(1:366,nF,'-m','linewidth',lw);
plot(find(nFLac~=0),nFLac(nFLac~=0),'xm');
plot(find(nFPreg~=0),nFPreg(nFPreg~=0),'om');
plot(1:366,nM,'-b','linewidth',lw);
ylabel('rec days');
xlabel('doy');
xlim([1 366]);
set(gca,'fontsize',14);
legend({'Total','Female','\rightarrowLactating','\rightarrowPregnant','Male'},'location','northwest');
legend box off;

% F/M histograms per-year
sexKey = {'F','M'};
ff(800,400);
for ii = 1:2
    subplot(1,2,ii);
    squirrels = strcmp(sqkey.sex,sexKey{ii});
    files = ~strcmp(sqkey.filename,'');
    years = sqkey.year(squirrels&files);
    n = histcounts(years);
    bar(unique(years),n);
    title(sprintf([sexKey{ii},', n = %i'],sum(n)));
end
