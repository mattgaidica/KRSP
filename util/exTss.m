function subTss = exTss(Tss,dts)

iRow = 1;
subTss = table;
while iRow <= numel(dts)
    curYear = year(dts(iRow));
    curMonth = month(dts(iRow));
    curDay = day(dts(iRow));
    % where else is this unique day?
    dtsIds = find(year(dts) == curYear & month(dts) == curMonth & day(dts) == curDay);
    tableRange = size(subTss,1)+1:size(subTss,1)+numel(dtsIds);
    subTss(tableRange,:) = repmat(Tss(Tss.year == curYear & Tss.month == curMonth & Tss.day == curDay,:),...
        [numel(tableRange),1]);
    iRow = dtsIds(end) + 1;
end