close all

filepath = '/Volumes/Seagate Expansion Drive 1/2019 AXY data/2019 Exported Axy';
files = dir(fullfile(filepath,'*.csv'));
% for iRow = 1:size(T,1)
for iRow = 1:numel(files)
%     filename = T.filename{iRow};
    filename = 'BJD15_21944_June-July2019_BJD15_2.csv';%files(iRow).name;
    saveFile = fullfile(filepath,[filename,'.jpg']);
    if isempty(filename) || isfile(saveFile)
        continue;
    end
    disp(filename);
    t = readtable(fullfile(filepath,filename));
    
    odba = abs(t.Var2)+abs(t.Var3)+abs(t.Var4);
    [Y,M,D] = datevec(T.DateStart(iRow));
    [h,m,s] = hms(T.TimeStart(iRow));
    dt_start = datetime(Y,M,D,h,m,s);
    [Y,M,D] = datevec(T.DateTurnedOff(iRow));
    [h,m,s] = hms(T.TimeTurnedOff(iRow));
    dt_end = datetime(Y,M,D,h,m,s);
    
    axy_dur = t.Var1(end) - t.Var1(1);
    human_dur = dt_end - dt_start;
    
    axy_diff = seconds(axy_dur - human_dur);
    
    h = ff(1200,400,2);
    plot(t.Var1,odba);
    hold on;
    plot(t.Var1,hour(t.Var1-hours(6)));
    ylabel('ODBA');
    yyaxis right;
    plot(t.Var1,t.Var5);
    ylabel('Collar temp (C)');
    hold on;
    plot([dt_start,dt_start],ylim,'g--','linewidth',2);
    plot([dt_end,dt_end],ylim,'r--','linewidth',2);
    xticks([datetime(year(t.Var1(1)),month(t.Var1(1)),day(t.Var1(1)),0,0,0)-days(1):...
        datetime(year(t.Var1(end)),month(t.Var1(end)),day(t.Var1(end)),0,0,0)+days(1)]);
    xlim auto;
    xtickangle(30);
    xlabel('date start (00:00)');
    set(gca,'fontsize',14);
    title({filename,sprintf('axy-human: %i seconds',axy_diff)},'interpreter','none');
    legend('ODBA','Hrs from 6AM','Collar Temp','Human Start','Human End');
    saveas(h,saveFile);
    close(h);
end