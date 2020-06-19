% load file
useFile = ''; % forces redo
dataPath = '/Users/matt/Documents/Data/KRSP/CompressedAxy';
if isempty(useFile)
    files = dir(fullfile(dataPath,'*.mat'));
    for iFile = 2:numel(files)
        saveFile = strrep(files(iFile).name,'.mat','_dn.txt');
        if ~isfile(fullfile(dataPath,saveFile)) % not analyzed
            useFile = files(iFile).name;
            break;
        end
    end
end
load(fullfile(dataPath,useFile));

warning('off','all');
dtDays = datetime(year(T.datetime),month(T.datetime),day(T.datetime));
allDates = unique(dtDays);
variable_names_types = [["filename", "string"]; ...
    ["day", "datetime"]; ...
    ["awake", "datetime"]; ...
    ["exit_nest", "datetime"]; ...
    ["enter_nest", "datetime"]; ...
    ["asleep", "datetime"]];
% Make table using fieldnames & value types from above
Tdn = table('Size',[numel(allDates),size(variable_names_types,1)],...
    'VariableNames', variable_names_types(:,1),...
    'VariableTypes', variable_names_types(:,2));
Tdn.day = allDates;
for iRow = 1:numel(allDates)
    Tdn.filename(iRow) = useFile;
end
% for iDay = 1:numel(allDates)
%     if isempty(Tmaster.day) || ~any(find(Tmaster.day == allDates(iDay)))
%         useRow = size(Tmaster,1)+1;
%         Tmaster.filename{useRow} = {useFile};
%         Tmaster.day(useRow) = allDates(iDay);
%     end
% end
warning('on','all');

rows = 2;
cols = 1;
clickColors = winter(4);
colors = lines(5);
op = 0.2;
lw = 1.5;
close all
clickStr = {'AWAKE','EXIT NEST','ENTER NEST','ASLEEP'};
for iDay = 1:3%numel(allDates)+1
    f = ff(1400,800);
    % data summary
    subplot(rows,cols,2);
    fmt = 'HH:mm';
    for jDay = 1:numel(allDates)
        dataIds_j = find(dtDays == allDates(jDay));
        t_j = linspace(T.datetime(dataIds_j(1)),T.datetime(dataIds_j(end)),24);
        yLow = secDay(T.datetime(dataIds_j(1)));
        yHigh = secDay(T.datetime(dataIds_j(end)));
        plot([jDay,jDay],[yLow,yHigh],'k:');
        hold on;
        ms = 15;
        plot(jDay,secDay(Tdn.awake(jDay)),'o','color',colors(5,:),'markersize',ms);
        plot(jDay,secDay(Tdn.exit_nest(jDay)),'x','color',colors(2,:),'markersize',ms);
        plot(jDay,secDay(Tdn.enter_nest(jDay)),'x','color',colors(2,:),'markersize',ms);
        plot(jDay,secDay(Tdn.asleep(jDay)),'o','color',colors(5,:),'markersize',ms);
    end
    xticks(1:numel(allDates));
    xticklabels(datestr(Tdn.day));
    xtickangle(30);
    ytickVals = linspace(1,86400,12);
    yticks(ytickVals);
    yticklabels(compose('%1.2f',ytickVals/3600));
    ylabel('hours');
    
    % day data
    subplot(rows,cols,1);
    if iDay <= numel(allDates)
        [sunrise,sunset,day_length] = sunriseSunset(allDates(iDay));
        dataIds = find(dtDays == allDates(iDay));
        t = linspace(T.datetime(dataIds(1)),T.datetime(dataIds(end)),24);
        
        yyaxis right;
        plot(T.datetime(dataIds),T.temp(dataIds),'color',colors(2,:),'linewidth',lw);
        set(gca,'ycolor',colors(2,:));
        ylabel('temp (C)');
        
        yyaxis left;
        plot(T.datetime(dataIds),T.odba(dataIds),'color',colors(5,:),'linewidth',lw);
        set(gca,'ycolor',colors(5,:));
        hold on;
        plot([sunrise,sunrise],[min(ylim),max(ylim)],'-','linewidth',5,'color',[colors(3,:) op]);
        plot([sunset,sunset],[min(ylim),max(ylim)],'-','linewidth',5,'color',[0 0 0 op]);
        xticks(t);
        xtickangle(30);
        ylabel('ODBA (g)');
        
        % click handler
        for iClick = 1:4
            title(clickStr{iClick});
            [x,y] = ginput(1);
            if x > 0
                clickDate = T.datetime(dataIds(1) + round(x * numel(dataIds)));
                plot([clickDate,clickDate],[min(ylim),max(ylim)],':','linewidth',2,'color',clickColors(iClick,:));
            else
                clickDate = NaN;
            end
            switch iClick
                case 1
                    Tdn.awake(iDay) = clickDate;
                case 2
                    Tdn.exit_nest(iDay) = clickDate;
                case 3
                    Tdn.enter_nest(iDay) = clickDate;
                case 4
                    Tdn.asleep(iDay) = clickDate;
            end
        end
    end
    
    if iDay > numel(allDates)
        answer = questdlg('Save results?','Day Night Analysis','Yes','No','No');
        if strcmp(answer,'Yes')
            writetable(Tdn,fullfile(dataPath,saveFile));
        end
    end
    
    close(f);
end