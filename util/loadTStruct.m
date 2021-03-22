function T = loadTStruct(iSq,sqkey,Tss)
loadPath = '/Users/matt/Box Sync/KRSP Axy Data/Temp';
if ~isnumeric(iSq) % iSq
    iSq = find(strcmp(iSq,sqkey.filename)); % overwrite
end
if ~isempty(iSq) && ~isempty(sqkey.filename{iSq}) && sqkey.isValid(iSq)
    load(fullfile(loadPath,sqkey.filename{iSq}));
    fprintf("%i/%i - %s\n",iSq,size(sqkey,1),sqkey.filename{iSq});
    T.datetime = T.datetime + minutes(sqkey.shiftMin(iSq));
    T.datetime.TimeZone = 'America/Whitehorse'; % automagically adds DST

    subTss = exTss(Tss,T.datetime);
    dls = subTss.day_length / 60; % min
    T = detect_sleepWake2(T,dls);
else
    fprintf("SKIPPING %i/%i - %s\n",iSq,size(sqkey,1),sqkey.filename{iSq});
    T = [];
end