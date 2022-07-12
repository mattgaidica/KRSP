if do
    sqkey = readtable('sqkey');
    filePath = '/Users/matt/Dropbox (University of Michigan)/from_box/KRSP Axy Data/Temp';
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
colors = lines(5); % females = 4, males = 5
rows = 3;
cols = 2;
ff(600,700,2);
subplot(rows,cols,[1,2]);
ms = 5;
plot(F_scatt_x,F_scatt_y,'.','markersize',ms,'color',colors(4,:));
hold on;
plot(M_scatt_x,M_scatt_y,'.','markersize',ms,'color',colors(5,:));
ylim([1 sqCount]);
xlim([1 366]);
xlabel('day of year');
ylabel('rec. session');
set(gca,'fontsize',14);
title(sprintf('%i recording sessions, %i recording days',sqCount,numel(M_scatt_x)+numel(F_scatt_y)));
grid on;
lns = [];
lns(1) = plot(-100,-100,'.','markersize',25,'color',colors(4,:));
lns(2) = plot(-100,-100,'.','markersize',25,'color',colors(5,:));
legend(lns,'Female','Male','location','southeast','box','off');

lw = 2;
subplot(rows,cols,[3,4]);
binEdges = linspace(0.5,366.5,366+1);
nM = histcounts(M_scatt_x,binEdges);
nF = histcounts(F_scatt_x,binEdges);
nFLac = histcounts(F_scattLac_x,binEdges);
nFPreg = histcounts(F_scattPreg_x,binEdges);
plot(1:366,nM+nF,'-k','linewidth',lw);
hold on;
plot(1:366,nF,'-','linewidth',lw,'color',colors(4,:));
plot(find(nFLac~=0),nFLac(nFLac~=0),'x','color',colors(4,:));
plot(find(nFPreg~=0),nFPreg(nFPreg~=0),'o','color',colors(4,:));
plot(1:366,nM,'-','linewidth',lw,'color',colors(5,:));
ylabel('rec. days');
xlabel('day of year');
xlim([1 366]);
set(gca,'fontsize',14);
legend({'Total','Female','\rightarrowLactating','\rightarrowPregnant','Male'},'location','northwest');
legend box off;
grid on;

% F/M histograms per-year
sexKey = {'F','M'};
sexLabels = {'Female','Male'};
for ii = 1:2
    subplot(rows,cols,ii+4);
    squirrels = strcmp(sqkey.sex,sexKey{ii});
    files = ~strcmp(sqkey.filename,'');
    years = sqkey.year(squirrels&files);
    n = histcounts(years,2013.5:1:2020.5);
    bar(2014:2020,n,'k');
    title(sprintf([sexLabels{ii},' (n = %i)'],sum(n)));
    ylim([0 70]);
    set(gca,'fontsize',14);
    xtickangle(30);
    grid on;
    ylabel('recording sessions');
end
saveas(gcf,fullfile(exportPath,'KRSP_axyDataCoverage.jpg'),'jpg');