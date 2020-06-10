sqkey = readtable('/Users/matt/Documents/Data/KRSP/sqkey.txt');
dataPath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';

colors = lines(4);
cs = [];
xs = [];
ys = [];
uniSq = 0;

startDoy = [];
useSquirrels = [];
for ii = 1:size(sqkey,1)
    if ~isempty(sqkey.filenames{ii})
        uniSq = uniSq + 1;
        useSquirrels(uniSq) = ii;
        fn = sqkey.filenames{ii};
        load(fullfile(dataPath,[fn,'.meta.mat']));
        startDoy(uniSq) = day(T_a.sunrise(1),'dayofyear');
        if size(T_a,1) > 1
            if T_a.sunrise(2) - T_a.sunrise(1) > day(2)
                disp(sqkey.filenames{ii})
            end
        end
    end
end

[~,k] = sort(startDoy);
useSquirrels = useSquirrels(k);

for ii = 1:numel(useSquirrels)
    fn = sqkey.filenames{useSquirrels(ii)};
    load(fullfile(dataPath,[fn,'.meta.mat']));
    for iDay = 1:size(T_a,1)
        xs = [xs;day(T_a.sunrise(iDay),'dayofyear')];
        ys = [ys;ii];
        if strcmp(sqkey.Sex(ii),'M')
            cs(numel(ys),:) = colors(1,:);
        elseif strcmp(sqkey.Sex(ii),'F')
            cs(numel(ys),:) = colors(2,:);
        else
            cs(numel(ys),:) = colors(4,:);
        end
    end
end

close all
lns = [];
ff;
scatter(xs,ys,20,cs,'filled');
hold on;
lns(1) = plot([0,0],[0,0],'.','color',colors(1,:),'markersize',20);
lns(2) = plot([0,0],[0,0],'.','color',colors(2,:),'markersize',20);
lns(3) = plot([0,0],[0,0],'.','color',colors(4,:),'markersize',20);
legend(lns,{'M','F','F-Lac'});
xlim([1 366]);
xlabel('doy');
ylim([1 uniSq]);
ylabel('squirrels');
yticks(ylim);

set(gca,'fontsize',14);