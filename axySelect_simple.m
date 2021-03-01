% !! set iSq manually
loadPath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
% sqkey = readtable('sqkey.txt');

Tss = readtable('/Users/matt/Documents/Data/KRSP/SunriseSunset/ss_2016.txt');
Tss_doys = day(Tss.sunrise,'dayofyear');
close all
while(iSq <= size(sqkey,1))
    if ~isempty(sqkey.filename{iSq})
        fprintf("Working on %s\n",sqkey.filename{iSq});
        load(fullfile(loadPath,sqkey.filename{iSq}));
    else
        fprintf("Skipping %s\n",sqkey.filename{iSq});
        sqkey.isValid(iSq) = false;
        iSq = iSq + 1;
        continue;
    end
    
    dayStart = find(round(secDay(T.datetime)/60) == 0,1,'first');
    theseOdba = [];
    theseOdba = T.odba(dayStart:end);
    theseOdba = [theseOdba;NaN(1440 - (numel(theseOdba) - 1440*(floor(numel(theseOdba)/1440))),1)];
    theseOdba = reshape(theseOdba,1440,[]);
    
    h = ff(1200,500);
    plot(theseOdba,'color',repmat(0.3,[4,1]));
    hold on;
    plot(nanmean(theseOdba,2)*4,'k-','linewidth',2); % scale this, easier to read
    plot(secDay(Tss.sunrise(day(T.datetime(1),'dayofyear')))/60,0,'go','markersize',20,'linewidth',2);
    plot(secDay(Tss.sunset(day(T.datetime(1),'dayofyear')))/60,0,'rx','markersize',20,'linewidth',2);
    xlim([1,1440]);
    title(sprintf("%03d: %s - %i days",iSq,sqkey.filename{iSq},size(theseOdba,2)),'interpreter','none');
    set(gca,'fontsize',16);
    
    text(10,max(ylim)/2,"Bad");
    text(1430,max(ylim)/2,"Good");
    
    [x,~] = ginput(1);
    if x > 0 && x < 720 % click left
        sqkey.isValid(iSq) = false;
    elseif x > 720 && x < 1440 % click right
        sqkey.isValid(iSq) = true;
    else
        close(h);
        break; % click margins
    end
    close(h);
    
    iSq = iSq + 1;
end

% these were compiled from the script below, need to fix/review them, but
% exclude for now
problemData = [84,85,88,179,181,190,194,197,199,247,294,348,357,366,417,418,534];
sqkey.isValid(problemData) = false;

% save sqkey!
writetable(sqkey,'sqkey.txt');
writetable(sqkey,'sqkey_isValid.txt');

%% check if day [sun ODBA] < [dark ODBA]
if false
    close all
    for iSq = 1:size(sqkey,1)
        if ~isempty(sqkey.filename{iSq}) & sqkey.filename{iSq}
            fprintf("Working on %s\n",sqkey.filename{iSq});
            load(fullfile(loadPath,sqkey.filename{iSq}));

            dayStart = find(round(secDay(T.datetime)/60) == 0,1,'first');
            theseOdba = [];
            theseOdba = T.odba(dayStart:end);
            theseOdba = [theseOdba;NaN(1440 - (numel(theseOdba) - 1440*(floor(numel(theseOdba)/1440))),1)];
            theseOdba = reshape(theseOdba,1440,[]);

            thisSunrise = round(secDay(Tss.sunrise(day(T.datetime(1),'dayofyear')))/60);
            thisSunset = round(secDay(Tss.sunset(day(T.datetime(1),'dayofyear')))/60);
            sunVals = theseOdba(thisSunrise:thisSunset,:);
            darkVals = [theseOdba(1:thisSunrise,:);theseOdba(thisSunset:end,:)];
            if nanmean(sunVals(:)) < nanmean(darkVals(:)) % suspect
                ff(1200,500);
                plot(theseOdba,'color',repmat(0.3,[4,1]));
                hold on;
                plot(nanmean(theseOdba,2)*4,'k-','linewidth',2); % scale this, easier to read
                plot(secDay(Tss.sunrise(day(T.datetime(1),'dayofyear')))/60,0,'go','markersize',20,'linewidth',2);
                plot(secDay(Tss.sunset(day(T.datetime(1),'dayofyear')))/60,0,'rx','markersize',20,'linewidth',2);
                xlim([1,1440]);
                title(sprintf("%03d: %s - %i days",iSq,sqkey.filename{iSq},size(theseOdba,2)),'interpreter','none');
                set(gca,'fontsize',16);
            end
        end
    end
end