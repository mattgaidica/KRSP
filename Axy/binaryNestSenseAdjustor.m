load('algorithmApp_session.mat');
sampleWindow = 60 * 30; % 60 * minutes
incBy = 10; % seconds

startFresh = true;
if exist('lastLoc','var')
    if lastLoc > 0
        answer = questdlg(sprintf("lastLoc = %i (%1.1f%%), reset?",lastLoc,100*lastLoc/numel(binNestSense)),...
            'Reset','Yes','No','Yes');
        if strcmp(answer,"No")
            startFresh = false;
        end
    end
end

if startFresh
    adj_binNestSense = binNestSense;
    lastLoc = 0;
end

colors = lines(5);
close all
ff(800,400);
while(1)
    transLocs = find(abs(diff(adj_binNestSense))==1);
    ii = find(transLocs > lastLoc,1,'first');
    if isempty(ii)
        break;
    end
    
    useRange = max([1,transLocs(ii)-sampleWindow+1]):min([numel(sense.smOdba),transLocs(ii)+sampleWindow]);
    
    yyaxis left;
    plot(useRange,odba(useRange),'-','color',colors(1,:),'linewidth',1);
    ylabel('ODBA');
    
    yyaxis right;
    plot(useRange,sense.smOdba(useRange),'k-','linewidth',3);
    hold on;
    plot(useRange,sense.tempZ(useRange),'r-','linewidth',3);
    binLine = plot(useRange,adj_binNestSense(useRange),'-','color',[colors(5,:),0.7],'linewidth',3);
    plot(useRange,binNestSense(useRange)/2+1.1,'-','color',[colors(5,:),0.15],'linewidth',3);
    set(gca,'ycolor','k');
    ylabel('Z');
    
    xlim([min(useRange),max(useRange)]);
    transSample = round((useRange(1)+useRange(end))/2);
    transLine = xline(transSample,'k:','linewidth',2);
    title(sprintf("Trans %i/%i",ii,numel(transLocs)));
    
    value = 0;
    ii_adj_binNestSense = adj_binNestSense;
    while(1) % enter, esc
        waitforbuttonpress;
        value = double(get(gcf,'CurrentCharacter'));
        if ismember(value,[13,27])
            break;
        end
        if value == 29 % right
            transSample = transSample + incBy;
        elseif value == 28 % left
            transSample = transSample - incBy;
        elseif value == 30 % up
            transSample = transLocs(ii); % reset
        end
        delete(transLine);
        transLine = xline(transSample,'k:','linewidth',2);
        
        adj_binNestSense(useRange) = ii_adj_binNestSense(useRange);
        if value == 31 % down
            % reject
            if ii < numel(transLocs)
                adj_binNestSense(transLocs(ii):transLocs(ii+1)) = ii_adj_binNestSense(transLocs(ii));
            end
            transSample = transLocs(ii); % reset
        else
            if transSample ~= transLocs(ii)
                if transSample > transLocs(ii)
                    adj_binNestSense(transLocs(ii):transSample) = ii_adj_binNestSense(transLocs(ii));
                else
                    adj_binNestSense(transSample:transLocs(ii)) = ~ii_adj_binNestSense(transLocs(ii));
                end
            end
        end
        delete(binLine);
        binLine = plot(useRange,adj_binNestSense(useRange),'-','color',[colors(5,:),0.7],'linewidth',3);
    end
    lastLoc = transSample;
    if value == 27
        break;
    end
    
    hold off;
    save("adjustSession.mat",'adj_binNestSense','lastLoc');
end
close all
%%
load("adjustSession.mat");
fprintf("Alikeness: %1.2f\n",sum(binNestSense==adj_binNestSense) / numel(adj_binNestSense));