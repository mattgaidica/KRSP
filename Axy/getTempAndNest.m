function [temp,nest] = getTempAndNest(temp,nSmooth)
temp = smoothdata(temp,'movmedian',30); % take care of dips
% find repeating values and fill them in
[C,~,ic] = unique(temp);
freq = accumarray(ic,1) / numel(temp);
repeatedValues = find(freq > 0.05);
for ii=1:numel(repeatedValues)
    fprintf("Removing repeated temperatures...\n");
    temp(temp==C(repeatedValues(ii))) = NaN;
end
temp = inpaint_nans(temp,4); % spring method
temp = smoothdata(temp,'gaussian',nSmooth);
nest = zeros(size(temp));
nDays = floor(numel(temp)/86400); % assumes 1Hz
fprintf("Performing k-means on temp:\n");
for iDay = 0:nDays-1
    fprintf("Day %i...\n",iDay);
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