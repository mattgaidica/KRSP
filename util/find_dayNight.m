% function find_dayNight()
% load('/Users/matt/Documents/Data/KRSP/ODBA/B6_Apr23_2015.csv.mat');
% load('/Users/matt/Documents/Data/KRSP/ODBA/B10_Apr23_2015.csv.mat');
% load('/Users/matt/Documents/Data/KRSP/ODBA/C19_Jan15_2017.csv.mat');
load('/Users/matt/Documents/Data/KRSP/ODBA/F5_Dec20_2015.csv.mat');
% load('/Users/matt/Documents/Data/KRSP/ODBA/D2_Apr23_2015.csv.mat');
% load('/Users/matt/Documents/Data/KRSP/ODBA/D17_Apr23_2015.csv.mat');
nSmooth = 360;
pThresh = 0.5;
colors = lines(3);

odba_smooth = smoothdata(T.odba,'gaussian',nSmooth);
odba_z = (odba_smooth - mean(odba_smooth)) / std(odba_smooth);
odba_bin = zeros(size(odba_z));
odba_bin(odba_z > 0) = 1;

all_awake_ids = find(diff(odba_bin) > 0);
all_asleep_ids = find(diff(odba_bin) < 0);

allDays = unique(datetime(year(T.datetime),month(T.datetime),day(T.datetime)));
sunrise_id = NaN(numel(allDays),1);
sunset_id = NaN(numel(allDays),1);
awake_id = NaN(numel(allDays),1);
asleep_id = NaN(numel(allDays),1);
for iDay = 1:numel(allDays)
    [sunrise,sunset] = sunriseSunset(allDays(iDay));
    tryId = find(T.datetime - sunrise > 0,1,'first');
    if ~isempty(tryId)
        sunrise_id(iDay) = tryId;
        [k,v] = closest(all_awake_ids,sunrise_id(iDay));
        awake_id(iDay) = v;
    end
    
    tryId = find(T.datetime - sunset > 0,1,'first');
    if ~isempty(tryId)
        sunset_id(iDay) = tryId;
        [k,v] = closest(all_asleep_ids,sunset_id(iDay));
        asleep_id(iDay) = v;
    end
end

op = 0.15;
fs = 12;
ly = 1.15;
lns = [];
close all
ff(1400,800);

yyaxis right;
hold on;
for iDay = 1:numel(allDays)
    if ~isempty(sunrise_id(iDay))
        lns(4) = plot([sunrise_id(iDay) sunrise_id(iDay)],ylim,'-','color',[colors(3,:)],'linewidth',4);
    end
    if ~isempty(sunset_id(iDay))
        lns(5) = plot([sunset_id(iDay) sunset_id(iDay)],ylim,'-','color',[0,0,0],'linewidth',4);
    end
end
lns(2) = plot(odba_z,'k','lineWidth',2);
hold on;
ylabel('z-odba');
set(gca,'ycolor','k')
xlabel_ids = sort([sunrise_id;sunset_id]);
xticks(xlabel_ids);
xticklabels(datestr(T.datetime(xlabel_ids)));
xtickangle(20);
xlabel('sunrise/sunset');
xlim([1 numel(odba_z)]);
set(gca,'fontsize',fs);
lns(3) = plot(odba_bin,'--r');

yyaxis left;
lns(1) = plot(T.odba,'color',[colors(1,:) op]);
ylabel('raw odba');
set(gca,'ycolor',colors(1,:))
set(gca,'fontsize',fs);

yyaxis right;
for iDay = 1:numel(allDays)
    if ~isempty(awake_id(iDay))
        text(awake_id(iDay),ly+0.05,{['#',num2str(iDay)],'\uparrow'},'fontsize',fs-2,'horizontalalign','center','color','r');
    end
    if ~isempty(asleep_id(iDay))
        text(asleep_id(iDay),ly,{['#',num2str(iDay)],'\downarrow'},'fontsize',fs-2,'horizontalalign','center','color','r');
    end
end

legend(lns,{'raw odba','gaussian z-odba','binary classifier','sunrise','sunset'},'location','eastoutside','fontsize',fs);

% % % % temp = smoothTemp(A.tempC(nRange),22.01);
% % % % temp = smoothdata(temp,'gaussian',nSmooth);
% % % % odba_raw = A.odba(nRange);
% % % % odba = smoothdata(odba_raw,'gaussian',nSmooth);
% % % % odba_z = (odba - mean(odba)) / std(odba);
% % % % 
% % % % s_odba = sort(odba);
% % % % odba_thresh = s_odba(round(numel(odba_raw)*pThresh));
% % % % 
% % % % close all
% % % % ff
% % % % 
% % % % subplot(221);
% % % % plot(A.tempC(nRange));
% % % % xlim([1 numel(nRange)]);
% % % % yyaxis right
% % % % plot(temp);
% % % % xlabel('temp C');
% % % % 
% % % % subplot(223);
% % % % plot(A.odba(nRange));
% % % % xlim([1 numel(nRange)]);
% % % % hold on;
% % % % yyaxis right
% % % % plot(odba);
% % % % xlabel('odba');
% % % % 
% % % % subplot(222);
% % % % odba_bin = zeros(size(odba));
% % % % odba_bin(odba_z>0) = 1;
% % % % plot(odba_z,'linewidth',2,'color','k');
% % % % hold on;
% % % % plot(odba_bin,'r');
% % % % xlim([1 numel(nRange)]);