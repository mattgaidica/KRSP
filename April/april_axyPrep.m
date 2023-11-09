function axy = april_axyPrep(axy,SN)
% Using movmedian for running median
axy.Xs = movmedian(axy.x, 91, 'omitnan');
axy.Ys = movmedian(axy.y, 91, 'omitnan');
axy.Zs = movmedian(axy.z, 91, 'omitnan');
axy.temp11 = movmedian(axy.temp, 11, 'omitnan');

% Grouping into 24-hour bins, assuming datetime is in datenum format
[~,~,axy.group] = histcounts(axy.datetime, 'BinWidth', days(1), 'BinLimits', [min(axy.datetime), max(axy.datetime)]);
% !! not sure what the code below does or protects
% % axy.group(1) = 1;
% % axy.group = fillmissing(axy.group, 'linear', 'SamplePoints', axy.datetime);

% Applying k-means clustering, assuming kmeansR is a user-defined function
G = findgroups(axy.group);
center = splitapply(@(x) april_kmeansR(x), axy.temp11, G);
center = fillmissing(center, 'linear');

% Inserting the centers back into the main array
for i = 1:length(center)
    axy.center(axy.group == i) = center(i);
end

% Calculating moving variance with a width of 420 and step of 30
axy.var = movvar(axy.temp, [419 0], 'omitnan'); % Note: Adjust the window as per requirement
axy.var(1) = 0;
try % !!not working with 10Hz data
    axy.var = fillmissing(axy.var, 'linear', 'SamplePoints', axy.datetime);
catch
    disp("Not filling missing var column in axyPrep");
end

% Assigning the SN variable
axy.Squirrel = repmat(SN, size(axy, 1), 1);
