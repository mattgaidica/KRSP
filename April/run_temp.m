warning('off', 'all');

newIds = [];
for ii = 276% 617:height(T_AxyDB)
    findId = find(contains(hashedAxyFiles.fullfile,T_AxyDB.md5_hash(ii)));
    if ~isempty(findId)
        fprintf("Loading row %i: %s\n",ii,T_AxyDB.md5_hash(ii));
        axy = readtable(hashedAxyFiles.fullfile(findId)); % ,'Delimiter',delimiter
        newHeaderLabels = strsplit(T_AxyDB.header_labels(ii),",");
        axy.Properties.VariableNames(1:numel(newHeaderLabels)) = newHeaderLabels;

        % handle case where date/time is string and combine them
        t_fmt = T_AxyDB.axy_t_fmt(ii);
        d_fmt = T_AxyDB.axy_d_fmt(ii);
        if ~any(contains(axy.Properties.VariableNames,"datetime"))
            if iscell(axy.date)
                if sum(axy.date{1}=='-') == 2
                    d_fmt = strrep(d_fmt,'/','-'); % replace
                end
                axy.date = datetime(axy.date,'Format',d_fmt);
            end
            if iscell(axy.time)
                axy.time = datetime(axy.time,'Format',t_fmt);
            end
            % combine date, time
            axy.datetime = axy.date + axy.time;
        end
        % !! fixing similar issue to above but requires datetime
        datetimeIncorrect = false;
        if iscell(axy.datetime)
            try
                if sum(axy.datetime{1}=='-') == 2
                    d_fmt = strrep(T_AxyDB.axy_d_fmt(ii),'/','-'); % replace
                end
                axy.datetime = datetime(axy.datetime,'Format',d_fmt+" "+t_fmt);
            catch
                datetimeIncorrect = true;
            end
        end
        if datetimeIncorrect
            disp("failed to convert datetime from string, skipping");
            continue;
        end
        axy.datetime.Format = 'dd-MMM-yyyy HH:mm:ss';

        if any(contains(axy.Properties.VariableNames,'datetime')) && any(contains(axy.Properties.VariableNames,'dtime'))
            if ~isdatetime(axy.datetime)
                newIds = [newIds;ii];
            else
                if seconds(max(diff(axy.datetime))) > 1
                    newIds = [newIds;ii];
                end
            end
        end
        % if (T_AxyDB.axy_fs == 10)
        %     newIds = [newIds;ii];
        % end
    end
end
newIds = unique(newIds);

warning('on', 'all');

%%

for ii = newIds
    disp(ii)
end