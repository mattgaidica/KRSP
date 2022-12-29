function [t,periOdba,periTemp] = periEventNest(binNestSense,odba,temp,nMin)

windowSamples = 60*nMin;
t = linspace(-nMin,nMin,windowSamples*2);

periOdba = {};
periTemp = {};
for ii = 1:2
    if ii==1
        useIds = find(diff(binNestSense) == 1);
    else
        useIds = find(diff(binNestSense) == -1);
    end
    kk = 0;
    periTemp{ii} = [];
    periOdba{ii} = [];
    for jj = 1:numel(useIds)
        useRange = useIds(jj)-windowSamples+1:useIds(jj)+windowSamples;
        if useRange(1) > 0 && useRange(end) < numel(binNestSense)
            kk = kk + 1;
            periTemp{ii}(kk,:) = temp(useRange); %#ok<*AGROW>
            periOdba{ii}(kk,:) = odba(useRange);
        end
    end
end