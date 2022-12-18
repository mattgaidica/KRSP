% function temp = filterTemp(T,nSmooth)
% find repeating values and fill them in
nSmooth = 120;
temp = T.temp;
[C,ia,ic] = unique(temp);
freq = accumarray(ic,1) / numel(temp);
repeatedValues = find(freq > 0.05); % 10%
for ii=1:numel(repeatedValues)
    temp(temp==C(repeatedValues(ii))) = NaN;
end
temp = inpaint_nans(temp,4); % spring method
temp = smoothdata(temp,'gaussian',nSmooth);
%%
nest = zeros(size(temp));
nDays = floor(numel(temp)/86400); % assumes 1Hz
for iDay = 0:nDays-1
    startIdx = 86400*iDay+1;
    if iDay == nDays
        endIdx = numel(temp);
    else
        endIdx = 86400*(iDay+1);
    end

    thisDay = temp(startIdx:endIdx);
    [IDX,C] = kmeans(thisDay,2);
    [~,nestIdx] = max(C);
    nest(find(IDX==nestIdx)+startIdx-1) = 1;
end
close all
ff(1200,400);
plot(nest);
hold on;
plot(strcmp(T.Nest,'Nest')+0.1,'g');
plot(normalize(T.odba,'range'),'k');
ylim([-1 2]);
yyaxis right;
plot(temp,'r-');