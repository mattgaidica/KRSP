
for ii = 1:height(T_AxyDB)
    % find hashed file
    findId = find(contains(hashedAxyFiles.fullfile,T_AxyDB.md5_hash(ii)));
    if ~isempty(findId)
        if any(T_AxyDB.deploy_date == NaT) || any(T_AxyDB.deploy_time == NaT)...
                || any(T_AxyDB.removed_date == NaT) || any(T_AxyDB.removed_time == NaT)
            disp(ii);
        end
    end
end