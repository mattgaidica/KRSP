function tempSorter_init(temp,nest,nMin,doAuto)
    windowHalfSize = 60*nMin;
    
    inNestMean = mean(temp(nest==1));
    inNestStd = std(temp(nest==1));
    outNestMean = mean(temp(nest==0));
    outNestStd = std(temp(nest==0));
    
    diff_nest = diff(nest);
    enterLocs = find(diff_nest==1);
    exitLocs = find(diff_nest==-1);
    enterNestTemp = NaN(1,windowHalfSize*2);
    exitNestTemp = NaN(1,windowHalfSize*2);

    jj = 0;
    for ii = 1:numel(enterLocs)
        useRange = enterLocs(ii)-windowHalfSize:enterLocs(ii)+windowHalfSize-1;
        if useRange(1) > 0 && useRange(end) < numel(temp)
            jj = jj + 1;
            enterNestTemp(jj,:) = temp(useRange);
        end
    end
    jj = 0;
    for ii = 1:numel(exitLocs)
        useRange = exitLocs(ii)-windowHalfSize:exitLocs(ii)+windowHalfSize-1;
        if useRange(1) > 0 && useRange(end) < numel(temp)
            jj = jj + 1;
            exitNestTemp(jj,:) = temp(useRange);
        end
    end
    t = linspace(-nMin,nMin,size(enterNestTemp,2));
    
    save(fullfile(pwd,'tempSession'),'t','enterLocs','exitLocs','enterNestTemp','exitNestTemp',...
        'inNestMean','inNestStd','outNestMean','outNestStd','doAuto');
    plotTrans(true);
end

function plotTrans(doInit)
    load(fullfile(pwd,'tempSession'),'t','enterNestTemp','exitNestTemp',...
        'inNestMean','inNestStd','outNestMean','outNestStd','doAuto');
    rows = 2;
    cols = 2;
    ylims = [min([enterNestTemp(:);exitNestTemp(:)]) max([enterNestTemp(:);exitNestTemp(:)])];
    ylabels = 'Temp (C)';
    xlabels = 'Time (min)';
    titleLabels = {'Enter Nest','Exit Nest'};
    fs = 14;
    colors = lines(max([size(enterNestTemp,1),size(exitNestTemp,1)]));
    if doInit
        close all;
        ff(1000,800);
        
        subplot(rows,cols,1);
        plot(t,enterNestTemp');
        xlim([min(t) max(t)]);
        ylim(ylims);
        ylabel(ylabels);
        set(gca,'fontsize',fs);
        title(titleLabels{1});
        grid on;
        
        if doAuto
            enter1 = drawrectangle('Label','ent_inc1','Color',[0 1 0],'Position',...
                [-20,min(ylim),20,(outNestMean+outNestStd)-min(ylim)]);
            enter2 = drawrectangle('Label','ent_inc2','Color',[0 1 0],'Position',...
                [0,inNestMean-inNestStd,20,max(ylim)-(inNestMean-inNestStd)]);
            enter3 = drawrectangle('Label','ent_exc1','Color',[1 0 0],'Position',...
                [-25,inNestMean-inNestStd*3,20,max(ylim)-(inNestMean-inNestStd*3)]);
            enter4 = drawrectangle('Label','ent_exc2','Color',[1 0 0],'Position',...
                [5,min(ylim),20,(outNestMean+outNestStd*3)-min(ylim)]);
        else
            diffx = mean(diff(xticks));
            diffy = mean(diff(yticks));
            enter1 = drawrectangle('Label','ent_inc1','Color',[0 1 0],'Position',...
                [min(xlim),min(ylim),diffx,diffy]);
            enter2 = drawrectangle('Label','ent_inc2','Color',[0 1 0],'Position',...
                [0,max(ylim)-diffy,diffx,diffy]);
            enter3 = drawrectangle('Label','ent_exc1','Color',[1 0 0],'Position',...
                [min(xlim),max(ylim)-diffy,diffx,diffy]);
            enter4 = drawrectangle('Label','ent_exc2','Color',[1 0 0],'Position',...
                [0,min(ylim),diffx,diffy]);
            addlistener(enter1,'ROIMoved',@allevents);
            addlistener(enter2,'ROIMoved',@allevents);
            addlistener(enter3,'ROIMoved',@allevents);
            addlistener(enter4,'ROIMoved',@allevents);
        end

        subplot(rows,cols,2);
        plot(t,exitNestTemp');
        xlim([min(t) max(t)]);
        ylim(ylims);
        ylabel(ylabels);
        set(gca,'fontsize',fs);
        title(titleLabels{2});
        grid on;
        
        if doAuto
            enter1 = drawrectangle('Label','exi_inc1','Color',[0 1 0],'Position',...
                [0,min(ylim),20,(outNestMean+outNestStd)-min(ylim)]);
            enter2 = drawrectangle('Label','exi_inc2','Color',[0 1 0],'Position',...
                [-20,inNestMean-inNestStd,20,max(ylim)-(inNestMean-inNestStd)]);
            enter3 = drawrectangle('Label','exi_exc1','Color',[1 0 0],'Position',...
                [5,inNestMean-inNestStd*3,20,max(ylim)-(inNestMean-inNestStd*3)]);
            enter4 = drawrectangle('Label','exi_exc2','Color',[1 0 0],'Position',...
                [-25,min(ylim),20,(outNestMean+outNestStd*3)-min(ylim)]);
        else
            diffx = mean(diff(xticks));
            diffy = mean(diff(yticks));
            exit1 = drawrectangle('Label','exi_inc1','Color',[0 1 0],'Position',...
                [min(xlim),max(ylim)-diffy,diffx,diffy]);
            exit2 = drawrectangle('Label','exi_inc2','Color',[0 1 0],'Position',...
                [0,min(ylim),diffx,diffy]);
            exit3 = drawrectangle('Label','exi_exc1','Color',[1 0 0],'Position',...
                [min(xlim),min(ylim),diffx,diffy]);
            exit4 = drawrectangle('Label','exi_exc2','Color',[1 0 0],'Position',...
                [0,max(ylim)-diffy,diffx,diffy]);
            addlistener(exit1,'ROIMoved',@allevents);
            addlistener(exit2,'ROIMoved',@allevents);
            addlistener(exit3,'ROIMoved',@allevents);
            addlistener(exit4,'ROIMoved',@allevents);
        end
        
        roiSave();
    else
        load(fullfile(pwd,'rois'),'roiTable');
        load(fullfile(pwd,'tempSession'),'t');
        exitEntryIdTable = table;
        
        subplot(rows,cols,3);
        inc1 = images.roi.Rectangle('Position',roiTable.pos{strcmp(roiTable.roi,"ent_inc1")});
        inc2 = images.roi.Rectangle('Position',roiTable.pos{strcmp(roiTable.roi,"ent_inc2")});
        exc1 = images.roi.Rectangle('Position',roiTable.pos{strcmp(roiTable.roi,"ent_exc1")});
        exc2 = images.roi.Rectangle('Position',roiTable.pos{strcmp(roiTable.roi,"ent_exc2")});
        useCount = 0;
        useIds = [];
        for ii = 1:size(enterNestTemp,1)
            theseData = enterNestTemp(ii,:);
            if (any(inROI(inc1,t,theseData)) && any(inROI(inc2,t,theseData))) && ...
                    ~(any(inROI(exc1,t,theseData)) || any(inROI(exc2,t,theseData)))
                plot(t,theseData,'-','color',colors(ii,:));
                hold on;
                useCount = useCount + 1;
                useIds(useCount) = ii;
            end
        end
        if ~isempty(useIds)
            ln1 = plot(t,mean(enterNestTemp(useIds,:)),'k','linewidth',2);
%             legend(ln1,{'avg'},'location','northwest','AutoUpdate','off');
        else
            legend off;
        end
        hold off;
        xlim([min(t) max(t)]);
        ylim(ylims);
        ylabel(ylabels);
        set(gca,'fontsize',fs);
        title(sprintf("%s (%i/%i)",titleLabels{1},useCount,ii));
        grid on;
        exitEntryIdTable.entries(1) = {useIds};
        
        subplot(rows,cols,4);
        inc1 = images.roi.Rectangle('Position',roiTable.pos{strcmp(roiTable.roi,"exi_inc1")});
        inc2 = images.roi.Rectangle('Position',roiTable.pos{strcmp(roiTable.roi,"exi_inc2")});
        exc1 = images.roi.Rectangle('Position',roiTable.pos{strcmp(roiTable.roi,"exi_exc1")});
        exc2 = images.roi.Rectangle('Position',roiTable.pos{strcmp(roiTable.roi,"exi_exc2")});
        useCount = 0;
        useIds = [];
        for ii = 1:size(exitNestTemp,1)
            theseData = exitNestTemp(ii,:);
            if (any(inROI(inc1,t,theseData)) && any(inROI(inc2,t,theseData))) && ...
                    ~(any(inROI(exc1,t,theseData)) || any(inROI(exc2,t,theseData)))
                plot(t,theseData,'-','color',colors(ii,:));
                hold on;
                useCount = useCount + 1;
                useIds(useCount) = ii;
            end
        end
        if ~isempty(useIds)
            ln1 = plot(t,mean(exitNestTemp(useIds,:)),'k','linewidth',2);
%             legend(ln1,{'avg'},'location','northwest','AutoUpdate','off');
        else
            legend off;
        end
        hold off;
        xlim([min(t) max(t)]);
        ylim(ylims);
        ylabel(ylabels);
        set(gca,'fontsize',fs);
        title(sprintf("%s (%i/%i)",titleLabels{2},useCount,ii));
        grid on;
        exitEntryIdTable.exits(1) = {useIds};
        
        save('exitEntryIdTable','exitEntryIdTable');
    end
end

function roiSave()
    ax = findobj(gcf, 'type', 'axes');
    roiTable = table;
    warning('off');
    rowCount = 0;
    for ii = 1:numel(ax)
        rects = findobj(ax(ii),'type','images.roi.Rectangle');
        for jj = 1:numel(rects)
            rowCount = rowCount + 1;
            roiTable.roi(rowCount) = string(rects(jj).Label);
            roiTable.pos(rowCount) = {rects(jj).Position};
        end
    end
    warning('on');
    save(fullfile(pwd,'rois'),'roiTable');
    plotTrans(false); % plot bottom
end

function allevents(src,evt)
    evname = evt.EventName;
    switch(evname)
        case{'ROIMoved'}
            roiSave();
    end
    
end