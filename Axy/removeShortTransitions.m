function nest = removeShortTransitions(nest,minSamples)

% perform recursively
while(1)
    diffNest = diff(nest);
    transLocs = find(abs(diffNest)==1);
    for ii = 1:numel(transLocs)-1
        jj = 0;
        if transLocs(ii+1)-transLocs(ii) < minSamples
            useRange = transLocs(ii)+1:transLocs(ii+1);
            nest(useRange) = repmat(nest(transLocs(ii)),[1,numel(useRange)]);
            fprintf("Removing %is interval at transition %i (%is)...\n",numel(useRange),ii,transLocs(ii));
            jj = jj+1;
        end
    end
    if jj == 0
        break;
    end
end