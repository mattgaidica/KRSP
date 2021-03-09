mastTitles = {'Mast','nMast'};
years_mast = [2014,2019]; % 2014,2019
years_nmast = [2015,2016,2017,2018,2020]; % 2015,2016,2017, *need 2018,2020
mast_years = {years_mast;years_nmast};

Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');

sIds = round(linspace(1,366,5));
seasonDoys = circshift(1:366,57);
useDoys = {seasonDoys(sIds(1):sIds(2)),seasonDoys(sIds(2):sIds(3)),...
    seasonDoys(sIds(3):sIds(4)),seasonDoys(sIds(4):sIds(5))};
monthNames = {'Winter','Spring','Summer','Autumn'};

close all
ff(1200,900);
for iSeason = 2:4
    mastIds = ismember(sq_doys,useDoys{iSeason}) & ismember(sq_years,years_mast);
    nmastIds = ismember(sq_doys,useDoys{iSeason}) & ismember(sq_years,years_nmast);
    
    sunrise = round(secDay(Tss.sunrise(round(mean(useDoys{iSeason})))) / 60);
    sunset = round(secDay(Tss.sunset(round(mean(useDoys{iSeason})))) / 60);
    
    %     subplot(3,2,prc(2,[iSeason-1,1]));
    %     plot(cumsum(circshift(mean(theseOdba_mast),sunrise)),'r-');
    %     hold on;
    %     plot(cumsum(circshift(mean(theseOdba_nmast),sunrise)),'k-');
    %     legend({'ODBA Mast','ODBA nMast'});
    %     xlim([1 1440]);
    
    subplot(3,1,iSeason-1);
    theseSleep_mast = cumsum(circshift(sq_asleep(mastIds,:),sunset,2),2);
    theseSleep_nmast = cumsum(circshift(sq_asleep(nmastIds,:),sunset,2),2);
    all_ps = [];
    for ii = 1:1440
        theseNmast = theseSleep_nmast(:,ii);
        theseMast = theseSleep_mast(:,ii);
        group = [zeros(size(theseNmast));ones(size(theseMast))];
        all_ps(ii) = anova1([theseNmast;theseMast],group,'off');
    end
    
    all_ps = pval_adjust(all_ps,'bonferroni');
    
    plot(mean(theseSleep_mast),'r-','linewidth',2);
    hold on;
    plot(mean(theseSleep_nmast),'k-','linewidth',2);
    plot(find(all_ps<0.001),0,'r*');
    legend({'Sleep Mast','Sleep nMast'},'location','northwest');
    xlim([1 1440]);
    
    title(monthNames{iSeason});
end

%%
sumSleep_mast = NaN(366,1);
sumSleep_nmast = NaN(366,1);
for iDoy = 1:366
    mastIds = sq_doys==iDoy & ismember(sq_years,years_mast);
    nmastIds = sq_doys==iDoy & ismember(sq_years,years_nmast);
    theseSleep_mast = sq_asleep(mastIds,:);
    theseSleep_nmast = sq_asleep(nmastIds,:);
    
    sunset = round(secDay(Tss.sunset(iDoy)) / 60);
    
    if numel(theseSleep_nmast)>1
        sumSleep_nmast(iDoy) = mean(sum(circshift(theseSleep_nmast,sunset,2),2));
    end
    if numel(theseSleep_mast)>1
        sumSleep_mast(iDoy) = mean(sum(circshift(theseSleep_mast,sunset,2),2));
    end
end

ff(1000,500);
lns = [];
lns(1) = plot(smoothdata(sumSleep_nmast,'gaussian',70),'k');
hold on;
lns(2) = plot(smoothdata(sumSleep_mast,'gaussian',70),'r');
title('Sum sleep');

for ii = 1:4
    plot([useDoys{ii}(1),useDoys{ii}(1)],ylim,'k:');
end

legend(lns,{'Sleep nMast','Sleep Mast'},'location','northwest');

%%
close all
ff(1400,500);
for iYear = 2014:2020
    sumSleep = NaN(366,1);
    for iDoy = 1:366
        theseIds = sq_doys==iDoy & sq_years==iYear;
        theseSleep = sq_asleep(theseIds,:);
        
        sumSleep(iDoy) = mean(mean(theseSleep,1));
    end
    if ismember(iYear,[2014,2019])
        dataColor = 'r';
    else
        dataColor = 'k';
    end
    plot(sumSleep,'linewidth',2,'color',dataColor);
    hold on;
end
for ii = 1:4
    plot([useDoys{ii}(1),useDoys{ii}(1)],ylim,'k:');
end
legend(compose("%i",2014:2020));
title('Sum sleep');
xlabel('day of year');
ylabel('mean asleep');
