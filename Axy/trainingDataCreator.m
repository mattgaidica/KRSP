samplesPerSegment = 60*60*8; % hours
close all
ff(1400,800);
nSegments = ceil(numel(temp)/samplesPerSegment);
if do
    transArr = [];
    jj = 0;
    lastLine = NaN;
    nestState = NaN;
    forceQuit = false;
    fs = 14;
    for ii = 1:nSegments
        startSample = (ii-1) * samplesPerSegment + 1;
        useRange = startSample:startSample+samplesPerSegment-1;
        tempRange = normalize(temp(useRange));
        odbaRange = normalize(odba(useRange),'range',[min(tempRange),max(tempRange)]);
        nestRange = nest(useRange);

        colors = lines(5);
        ylabel('raw axy');
        plot(useRange,odbaRange,'k-','linewidth',1);
        hold on;
        plot(useRange,tempRange,'r-','linewidth',1);
    %     ylabel('ODBA');

        plot(useRange,nestRange,'-','color',[colors(5,:),0.35],'linewidth',2);
        hold on;
    %     set(gca,'ycolor',colors(5,:));
        set(gca,'fontsize',fs);
    %     ylim([-.1 1.1]);
        yticks([0 1]);
        yticklabels({'Out of Nest','In Nest'});
        ylim([-max(abs(ylim)) max(abs(ylim))]);
        xlim([min(useRange),max(useRange)]);
        xlabel('Samples [Click to Quit]');
        grid on;
        legend({'ODBA','Nest','Temp'},'AutoUpdate','off')
        title(sprintf("Segment %i/%i",ii,nSegments));
        if ii == 1
            text(min(xlim)-60,max(ylim),'Undo','horizontalalignment','right','fontsize',fs+2);
            text(max(xlim)+60,max(ylim),'Next','horizontalalignment','left','fontsize',fs+2);
        end
        lns = [];
        kk = 0;
        while(1)
            [x,y] = ginput(1);
            if x > useRange(end) % next
                break;
            end
            if x < useRange(1) % undo
                delete(lns(kk));
                if kk > 1
                    kk = kk - 1;
                end
            end
            if y < min(ylim) % quit
                forceQuit = true;
                break;
            end

            if isnan(nestState)
                if y > 0
                    nestState = 1;
                else
                    nestState = 0;
                end
            end
            if nestState == 0
                markerStyle = 'ro';
            else
                markerStyle = 'go';
            end
            jj = jj + 1;
            kk = kk + 1;
            lns(kk) = plot(round(x),nestState,markerStyle,'markersize',10);
            transArr(jj,:) = [round(x),nestState];
            nestState = ~nestState;
        end
        hold off;
        if forceQuit
            break;
        end
    end
end
%% make nest
trainingNest = NaN(size(nest));
nestState = ~transArr(1,2);
for ii = 1:numel(trainingNest)
    useIdx = find(transArr(:,1) == ii);
    if ~isempty(useIdx)
        nestState = transArr(useIdx,2);
    end
    trainingNest(ii) = nestState;
end

close all;
ff(1200,400);
plot(odba,'k-');
hold on;
plot(trainingNest,'g-');
plot(nest+1.1,'m-');

fprintf('%1.2f alike\n',sum(trainingNest(1:ii)==nest(1:ii))/ii);
