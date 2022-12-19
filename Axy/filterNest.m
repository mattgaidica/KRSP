function newNest = filterNest(nest,odba,temp)
% use unshiftedNest
tempGrad = gradient(temp);
odbaDir = normalize(odba);
nestSense = temp;
nestSense(sign(tempGrad)~=sign(-odbaDir)) = NaN;
binNestSense = zeros(size(nestSense));
binNestSense(~isnan(nestSense)) = 1;
probCorrect = smoothdata(binNestSense.*sign(tempGrad)+normalize(tempGrad,'range',[-1 1]),'gaussian',60*5);

probThresh = 0.5;
diffNest = diff(nest);
transEnter = find(diffNest == 1 & probCorrect(2:end) > probThresh);
transExit = find(diffNest == -1 & probCorrect(2:end) < -probThresh);
transIds = sort([1;transEnter;transExit;numel(nest)]);
newNest = NaN(size(nest));
nestState = ~nest(2);
lastValid = 1;
for ii = 2:numel(transIds)
    if nest(transIds(ii)) ~= nestState
        nestState = ~nestState;
        useRange = transIds(lastValid):transIds(ii);
        newNest(useRange) = repmat(nestState,[1,numel(useRange)]);
        lastValid = ii;
    end
end


nestNaN = nest;
nestNaN(abs(probCorrect)<probThresh) = NaN;
colors = lines(5);
close all
ff(1200,400);
plot(normalize(odba,'range')*3,'color',[0 0 0 0.4]);
hold on;
plot(nest,'-','color',colors(5,:),'linewidth',2);
ylim([-1 3]);
plot(probCorrect,'b-');
plot(nestNaN,'-','color','g','linewidth',3);
yyaxis right;
hold on;
plot(temp,'r:');
plot(nestSense,'r-','linewidth',2);
set(gca,'ycolor','r');
grid on;
yyaxis left;
plot(newNest+1.1,'-','color',[colors(3,:),.5],'linewidth',3);