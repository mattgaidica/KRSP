% fix NaN squirrel ids
if ~exist('squirrel','var')
    squirrel = readtable('/Users/matt/Documents/MATLAB/KRSP/R/krsp_squirrel_tagsOnly.csv','Format','%s%s%s');
    squirrel.squirrel_id = str2double(squirrel.id);
end
for iSq = 1:size(sqkey,1)
    if isnan(sqkey.squirrel_id(iSq)) && strcmp(sqkey.tag_left(iSq),'') == 0
        if strcmp(sqkey.tag_left(iSq),'-') == 1
            id = find(strcmp(sqkey.tag_right(iSq),squirrel.tagrt) & strcmp('NULL',squirrel.taglft));
        elseif strcmp(sqkey.tag_right(iSq),'-') == 1
            id = find(strcmp(sqkey.tag_left(iSq),squirrel.taglft) & strcmp('NULL',squirrel.tagrt));
        else
            id = find(strcmp(sqkey.tag_left(iSq),squirrel.taglft) & strcmp(sqkey.tag_right(iSq),squirrel.tagrt));
        end
        if numel(id) == 1
            disp(squirrel.squirrel_id(id));
            sqkey.squirrel_id(iSq) = squirrel.squirrel_id(id);
        end
    end
end

% manual fix to iSq
% where tag_right = {'A7240'}, enter tag_left = 'M5841'
% where tag_left = 'M2264', enter tag_right = 'M7804'

%%

% litter = readtable('/Users/matt/Documents/MATLAB/KRSP/R/krsp_litter.csv');
% juvenile = readtable('/Users/matt/Documents/MATLAB/KRSP/R/krsp_juvenile.csv');
% growth = readtable('/Users/matt/Documents/MATLAB/KRSP/R/krsp_growth.csv');
for iSq = 1:size(sqkey,1)
    if isempty(sqkey.filename{iSq}) || sqkey.isValid(iSq) == 0
        continue;
    end
    T = loadTStruct(iSq,sqkey,Tss);
    if isempty(T)
        continue;
    end
    iValid = iValid + 1;
    
    sqkey.rec_startDate(iSq) = T.datetime(1);
    sqkey.rec_endDate(iSq) = T.datetime(end);
    sqkey.rec_midDate(iSq) = T.datetime(1) + hours(T.datetime(end)-T.datetime(1)) / 2;
    sqkey.rec_midDoy(iSq) = day(sqkey.rec_midDate(iSq),'dayofyear');
    sqkey.rec_mins(iSq) = size(T,1);
end

%% preg status
for iSq = 1:size(sqkey,1)
    if any(strcmp(sqkey.sex_status{iSq},{'lactating'})) % 'Pre-pregnancy','pregnant'
        sqkey.is_preg(iSq) = 1;
    else
        sqkey.is_preg(iSq) = 0;
    end
end

%%
litter.fieldBDate.TimeZone = 'America/Whitehorse';
for iSq = 1:size(sqkey,1)
    sqkey.rec_litterId(iSq) = NaN;
    sqkey.rec_litterSize(iSq) = NaN;
    if sqkey.is_preg(iSq) == 1
        rowIds = find(litter.squirrel_id == sqkey.squirrel_id(iSq));
        if ~isempty(rowIds)
            % find nearest litter
            litterDaysFromRec = days(litter.fieldBDate(rowIds) - sqkey.rec_startDate(iSq));
            [v,k] = min(litterDaysFromRec);
            
            litterId = litter.id(rowIds(k));
            sqkey.rec_litterId(iSq) = litterId;
            juvs = find(juvenile.litter_id == litterId);
            sqkey.rec_litterSize(iSq) = numel(juvs);
        end
    end
end
% >> writetable(sqkey,'sqkey');