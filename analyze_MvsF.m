load('sqkey.mat');
dataPath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
useDoy = true;

colors = lines(4);
cs = [];
xs = [];
ys = [];
uniSq = 0;

sortDates = [];
useSquirrels = [];
for ii = 1:size(sqkey,1)
    if ~isempty(sqkey.filenames{ii})
        uniSq = uniSq + 1;
        useSquirrels(uniSq) = ii;
        fn = sqkey.filenames{ii};
        fn_meta = strrep(fn,'.mat','_meta.mat');
        load(fullfile(dataPath,fn_meta));
        if useDoy
            sortDates(uniSq) = day(T_a.sunrise(1),'dayofyear');
        else
            sortDates(uniSq) = datenum(T_a.sunrise(1));
        end
        if size(T_a,1) > 1 % BAD DATES
            if T_a.sunrise(2) - T_a.sunrise(1) > day(2)
                disp(sqkey.filenames{ii})
            end
        end
    end
end

[~,k] = sort(sortDates);
useSquirrels = useSquirrels(k);

warning('off','all');
T_scatt = table;
dayCount = 0;
for ii = 1:numel(useSquirrels)
    fn = sqkey.filenames{useSquirrels(ii)};
    fn_meta = strrep(fn,'.mat','_meta.mat');
    load(fullfile(dataPath,fn_meta));
    for iDay = 1:size(T_a,1)
        dayCount = dayCount + 1;
        if useDoy
            T_scatt.xs(dayCount) = day(T_a.sunrise(iDay),'dayofyear');
        else
            T_scatt.xs(dayCount) = T_a.sunrise(iDay);
        end
        T_scatt.ys(dayCount) = ii;
        if strcmp(sqkey.Sex(ii),'M')
            T_scatt.cs(dayCount) = 1;
        elseif strcmp(sqkey.Sex(ii),'F')
            T_scatt.cs(dayCount) = 2;
        else
            T_scatt.cs(dayCount) = 4;
        end
    end
end
warning('on','all');

close all
lns = [];
ff;
scatter(T_scatt.xs,T_scatt.ys,20,colors(T_scatt.cs,:),'filled');
hold on;
if useDoy
    lns(1) = plot([0,0],[0,0],'.','color',colors(1,:),'markersize',20);
    lns(2) = plot([0,0],[0,0],'.','color',colors(2,:),'markersize',20);
    lns(3) = plot([0,0],[0,0],'.','color',colors(4,:),'markersize',20);
    xlim([1 366]);
    xlabel('day of year');
else
    lns(1) = plot([T_scatt.xs(1),T_scatt.xs(1)],[0,0],'.','color',colors(1,:),'markersize',20);
    lns(2) = plot([T_scatt.xs(1),T_scatt.xs(1)],[0,0],'.','color',colors(2,:),'markersize',20);
    lns(3) = plot([T_scatt.xs(1),T_scatt.xs(1)],[0,0],'.','color',colors(4,:),'markersize',20);
    xlim([min(T_scatt.xs) max(T_scatt.xs)]);
end
legend(lns,{'M','F','F-Lac'},'location','best');
ylim([1 uniSq]);
ylabel('squirrels');
yticks(ylim);

set(gca,'fontsize',20);