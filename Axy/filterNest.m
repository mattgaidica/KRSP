function newNest = filterNest(nest,odba,temp)
doPlot = true;

probThresh = 0.5;
% use unshiftedNest
tempGrad = gradient(temp);
odbaDir = normalize(odba);
nestSense = temp;
nestSense(sign(tempGrad)~=sign(-odbaDir)) = NaN;
binNestSense = zeros(size(nestSense));
binNestSense(~isnan(nestSense)) = 1;
probCorrect = smoothdata(binNestSense.*sign(tempGrad)+normalize(tempGrad,'range',[-1 1]),'gaussian',60);
nestNaN = nest;
nestNaN(abs(probCorrect)<probThresh) = NaN;
newNest = round(inpaint_nans(nestNaN,5));

% below does not work because it misses trans out of nest too often
% % diffNest = diff(nest);
% % % only use probable transitions in the right direction (probability is +/-)
% % transEnter = find(diffNest == 1 & probCorrect(2:end) > probThresh);
% % transExit = find(diffNest == -1 & probCorrect(2:end) < -probThresh);
% % transIds = sort([1;transEnter;transExit;numel(nest)]);
% % newNest = NaN(size(nest));
% % nestState = ~nest(2);
% % lastValid = 1;
% % for ii = 2:numel(transIds)
% %     if nest(transIds(ii)) ~= nestState
% %         nestState = ~nestState;
% %         useRange = transIds(lastValid):transIds(ii);
% %         newNest(useRange) = repmat(nestState,[1,numel(useRange)]);
% %         lastValid = ii;
% %     else
% %         fprintf("Bad transition, filling in...\n");
% %     end
% % end

if doPlot
    colors = lines(5);
    close all
    ff(1200,400);
    plot(normalize(odba,'range')*3,'color',[0 0 0 0.4]);
    hold on;
    plot(nest,'-','color',colors(5,:),'linewidth',2);
    ylim([-2 3]);
    plot(probCorrect,'b-');
    plot(nestNaN,'-','color',[0 1 0 0.5],'linewidth',3);
    yyaxis right;
    hold on;
    plot(temp,'r:');
    plot(nestSense,'r-','linewidth',2);
    set(gca,'ycolor','r');
    grid on;
    yyaxis left;
    plot(newNest+1.1,'-','color',[colors(3,:),.5],'linewidth',3);
    xlim([0 86400/12]);
end