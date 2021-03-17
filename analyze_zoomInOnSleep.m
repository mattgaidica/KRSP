if do
    A = readtable('/Users/matt/Box Sync/KRSP Axy Data/Dantzer Lab Axy Data/2019 AXY data/Axy Data/BJD10/21961/MAY 2019/BJD10_1.csv');
    load('/Users/matt/Box Sync/KRSP Axy Data/Temp/BJD10_21961_MAY2019_BJD10_1_nest_trim_nest__20190505.mat');
    T = detect_sleepWake2(T,60);
    Tawake = make_Tawake(T);
end

%% look at brief awakenings to determine nature of behavior
% setup here
cumAwake = zeros(size(T,1),1);
awakeIdx = find(diff(T.awake) == 1)+1;
for iRow = awakeIdx'
    useIdx = iRow;
    curSum = 1;
    while (useIdx < size(T,1) & T.awake(useIdx) == 1)
        cumAwake(useIdx) = curSum;
        curSum = curSum + 1;
        useIdx = useIdx + 1;
    end
end

%% continuous measure on A
x_movMean = movmean(A.Var2,60);
x_grad = gradient(A.Var2);
x_cumsum = cumsum(x_grad) - x_movMean;
x_cumsum_abs = cumsum(abs(x_grad)) - x_movMean;

close all
ff(1200,600);
plot(A.Var2);
hold on;
plot(x_cumsum);
yyaxis right;
plot(x_cumsum);
%% simulate axy
t = linspace(0,2*pi,1000);
sig1 = sin(t);
sig2 = [zeros(1,500),ones(1,250),zeros(1,250)];
sig3 = [zeros(1,500),ones(1,500)*0.25];
allSigs = [sig1;sig2;sig3];

plot(sig1);
hold on;
plot(sig2);
plot(sig3);

rows = 3;
cols = 2;
close all
ff(1400,500);
for ii=1:3
    subplot(rows,cols,prc(cols,[ii,1]));
    x = allSigs(ii,:) - allSigs(ii,1);
    plot(x); hold on;
    cumx = cumsum(gradient(x));
    cumAbsx = cumsum(abs(gradient(x)));
    subplot(rows,cols,prc(cols,[ii,2]));
    plot(cumx); hold on;
    plot(cumAbsx);
end
%% 
% Say we are trying to identify REM awake/small position shifts: these CAN
% be associated with small, high-amplitude movements resulting in high
% ODBA. This behavior would allow changes in orientation (turning), but not
% large changes in distance. 
[locs,pks] = peakseek(cumAwake,1,1);
upToMinutes = 10;
viewExtraMinutes = 1;
if true
    shortLocs = locs(pks<=upToMinutes);
    shortPks = pks(pks<=upToMinutes);
else
    shortLocs = locs(pks>60)-10;
    shortPks = pks(pks>60)-10;
end
close all
clc
doPlot = false;
sigCount = 0;
for shortEnd = shortLocs % 2 is sig, 4 is clearly movement, 6 should be sig
    AIdx = find(A.Var1 == T.datetime(shortEnd));
    
    shortStart = find(cumAwake(1:shortEnd)==0,1,'last')+1;
    showTRange = (shortStart-viewExtraMinutes):(shortEnd+viewExtraMinutes);
    shortLength = shortEnd - shortStart;

    % handle Fs!
    analyzeRawRange = (AIdx-60*shortLength-1):AIdx;
    showRawRange = (AIdx-60*(shortLength+viewExtraMinutes)+1):(AIdx+viewExtraMinutes*60);
    xShow = A.Var2(showRawRange) - A.Var2(showRawRange(1));
    yShow = A.Var3(showRawRange) - A.Var3(showRawRange(1));
    zShow = A.Var4(showRawRange) - A.Var4(showRawRange(1));
    
    if do
        nSurr = 10000;
        surrScores = NaN(nSurr,1);
        locPool = find(cumAwake(1:end-shortLength)>upToMinutes);
        for iSurr = 1:nSurr
            randWakeLoc = randsample(locPool,1);
            surrStart = find(A.Var1 == T.datetime(randWakeLoc));
            surrRange = surrStart:surrStart+(upToMinutes*60)-1;
            x = A.Var2(surrRange) - A.Var2(surrRange(1));
            y = A.Var3(surrRange) - A.Var3(surrRange(1));
            z = A.Var4(surrRange) - A.Var4(surrRange(1));

            cumx = cumsum(gradient(x));
            cumy = cumsum(gradient(y));
            cumz = cumsum(gradient(z));

            cumAbsx = cumsum(abs(gradient(x)));
            cumAbsy = cumsum(abs(gradient(y)));
            cumAbsz = cumsum(abs(gradient(z)));

            surrDist = cumAbsx(end)+cumAbsy(end)+cumAbsz(end);
            surrOrient = abs(cumx(end))+abs(cumy(end))+abs(cumz(end));
            surrScores(iSurr) = (surrDist * surrOrient) / upToMinutes;

    % % % %         h = ff(700,900);
    % % % %         subplot(311);
    % % % %         plot(x); hold on
    % % % %         plot(cumx);
    % % % %         plot(cumAbsx);
    % % % %         subplot(312);
    % % % %         plot(A.Var2);
    % % % %         hold on;
    % % % %         xline(surrRange(1),'g-');
    % % % %         xline(surrRange(end),'r-');
    % % % %         xlim([surrRange(1)-60000 surrRange(end)+60000]);
    % % % %         title(sprintf('Score: %1.4f',surrScores(iSurr)));
    % % % %         subplot(313);
    % % % %         plot(A.Var2);
    % % % %         hold on;
    % % % %         xline(surrRange(1),'g-');
    % % % %         xline(surrRange(end),'r-');
    % % % %         xlim([surrRange(1)-1000 surrRange(end)+1000]);
    % % % %         title(sprintf('Score: %1.4f',surrScores(iSurr)));
    % % % %         close(h)
        end
        do = false;
    end
    
    x = A.Var2(analyzeRawRange) - A.Var2(analyzeRawRange(1));
    y = A.Var3(analyzeRawRange) - A.Var3(analyzeRawRange(1));
    z = A.Var4(analyzeRawRange) - A.Var4(analyzeRawRange(1));
    
    cumx = cumsum(gradient(x));
    cumy = cumsum(gradient(y));
    cumz = cumsum(gradient(z));
    
    cumAbsx = cumsum(abs(gradient(x)));
    cumAbsy = cumsum(abs(gradient(y)));
    cumAbsz = cumsum(abs(gradient(z)));
    
    thisDist = cumAbsx(end)+cumAbsy(end)+cumAbsz(end);
    thisOrient = abs(cumx(end))+abs(cumy(end))+abs(cumz(end));
    thisScore = (thisOrient * thisDist) / shortLength;
    p = 1 - sum(thisScore < surrScores)/nSurr;
    if p < 0.05 && p >= 0.01
        fprintf('*');
        sigCount = sigCount + 1;
    elseif p < 0.01
        fprintf('**');
        sigCount = sigCount + 1;
    end
    fprintf('p = %1.4f, score: %1.4f\n',thisScore,p)
    if doPlot
        rows = 2;
        cols = 8;
        ff(1400,500);
        subplot(rows,cols,[1,2]);
        plot(T.odba(showTRange));
        yyaxis right;
        plot(T.awake(showTRange),'k');
        set(gca,'ycolor','k');
        ylim([-0.5 1.5]);
        grid on;
        xlabel('Minutes (1-min res.)');
        xlim(size(showTRange));
        xticks(min(xlim):max(xlim));
        xticklabels(xticks-1); % 0-based
        title('T-struct');
        legend({'ODBA','Awake'},'autoupdate','off','location','westoutside');
        xline(min(xlim)+viewExtraMinutes,'k--');
        xline(max(xlim)-viewExtraMinutes,'k--');

        subplot(rows,cols,[9,10]);
        t = linspace(0,numel(showTRange),numel(xShow));
        plot(t,xShow);
        hold on;
        plot(t,yShow);
        plot(t,zShow);
        xlabel('Minutes (1-sec res.)');
        xlim(size(showTRange));
        xticks(min(xlim):max(xlim));
        xticklabels(xticks-1); % 0-based
        title('Raw Axy');
        legend({'X','Y','Z'},'autoupdate','off','location','westoutside');
        grid on;
        c = colorbar('location','southoutside');
        colormap(jet);
        c.Ticks = [0,1];
        c.TickLabels = {'Start','End'};
        xline(min(xlim)+viewExtraMinutes,'k--');
        xline(max(xlim)-viewExtraMinutes,'k--');

        colors = jet(numel(analyzeRawRange));
        subplot(rows,cols,[3,4,11,12]);
        plot3(cumx,cumy,cumz,'color',repmat(0.7,[1 3]),'linewidth',1.5);
        hold on;
        scatter3(cumx,cumy,cumz,35,colors,'filled');
        plot3([cumx(1),cumx(end)],[cumy(1),cumy(end)],[cumz(1),cumz(end)],'r','linewidth',3);
        xlabel('X');
        xlabel('Y');
        xlabel('Z');
        grid on;
        title(sprintf('Orient dist \nd = %2.2f',thisOrient));

        subplot(rows,cols,[5,6,13,14]);
        plot3(cumAbsx,cumAbsy,cumAbsz,'color',repmat(0.7,[1 3]),'linewidth',1.5);
        hold on;
        scatter3(cumAbsx,cumAbsy,cumAbsz,35,colors,'filled');
        plot3([cumAbsx(1),cumAbsx(end)],[cumAbsy(1),cumAbsy(end)],[cumAbsz(1),cumAbsz(end)],'r','linewidth',3);
        xlabel('X');
        xlabel('Y');
        xlabel('Z');
        grid on;
        
        title(sprintf('Move dist \nd = %2.2f',thisDist));

        subplot(rows,cols,[7,8,15,16]);
        histogram(surrScores,1:100);
        hold on;
        xline(thisScore,'r-','linewidth',2);

        title(sprintf('Score: %1.4f\nLess than most? p = %1.4f',thisScore,p));
        xlabel('All Scores');
        ylabel('count');
    end
end
fprintf('\n%i/%i p < 0.05 — %2.1f%% of <%i min awake bouts are actually sleep\n',sigCount,numel(shortLocs),100*sigCount/numel(shortLocs),upToMinutes);

%% apply algorithm to entire raw table
Fs = 1; % 1 Hz
windowSize = Fs * 60;
distOrient = NaN(size(A,1),1); % should I inpaint the ends or handle NaN later?
distMoved = NaN(size(A,1),1);
for ii = 0.5*10e5:0.7*10e5%halfWindow:size(A,1)-halfWindow
    nRange = (ii-windowSize+1):ii;
    x = A.Var2(nRange) - A.Var2(ii);
    y = A.Var3(nRange) - A.Var3(ii);
    z = A.Var4(nRange) - A.Var4(ii);
    
    cumx = cumsum(gradient(x));
    cumy = cumsum(gradient(y));
    cumz = cumsum(gradient(z));
    
    cumAbsx = cumsum(abs(gradient(x)));
    cumAbsy = cumsum(abs(gradient(y)));
    cumAbsz = cumsum(abs(gradient(z)));
    
    distOrient(ii) = abs(cumx(end)) + abs(cumy(end)) + abs(cumz(end));
    distMoved(ii) = cumAbsx(end) + cumAbsy(end) + cumAbsz(end);
end

%% show raw data
nSmooth = 91;
X = A.Var2;
Y = A.Var3;
Z = A.Var4;
Aodba = abs(X-medfilt1(X,nSmooth))...
    + abs(Y-medfilt1(Y,nSmooth))...
    + abs(Z-medfilt1(Z,nSmooth));

ff(1200,800);
plot(Aodba,'g-');
yyaxis right;
plot(normalize(distMoved,'range'),'b-');
hold on;
plot(normalize(distOrient,'range'),'r-');
plot(normalize(distMoved.*distOrient,'range'),'k-','linewidth',2);
legend({'ODBA','moved','orient','score'});
xlim([5.0005e+05 5.0005e+05 + 60*5]);